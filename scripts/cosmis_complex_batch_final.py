import csv
import gzip
import json
import re
import sys
# Adding directory to python path
# Needs a better fix 
sys.path.append('/wynton/home/rotation/gyasu/capra_lab/cosmis')

import logging
import numpy as np
from argparse import ArgumentParser
from collections import defaultdict
import os

from Bio import SeqIO
from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils import seq1
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO


from cosmis import seq_utils, pdb_utils

# Terminal line command
'''
python scripts/cosmis_complex_batch_final.py -c data_paths.json -i /wynton/home/rotation/gyasu/FoldDock/output_txt_dir/ -o AFcomplex_chainA.tsv -l AFcomplex_chainA.log --chain A -d /wynton/home/rotation/gyasu/FoldDock/output_pdbs/

'''
def parse_cmd():
    """
    Specifies command-line flags and parses command-line options.

    Returns
    -------
    ArgumentParser
        An object of type ArgumentParser containing the parsed command-line
        arguments.

    """
    parser = ArgumentParser()
    parser.add_argument('-c', '--config', dest='config', required=True,
                        type=str, help='A JSON file specifying options.')
    parser.add_argument('-i', '--input', dest='input', required=True,
                        type=str, help='''Directory containing text files where each file 
                        contains pairs of PDB chain ID:Ensembl transcript ID associated 
                        with the input PDB file.''')
    parser.add_argument('-o', '--output', dest='output_file', required=True,
                        type=str, help='''Output file to store the COSMIS scores
                        of the protein.''')
    parser.add_argument('-l', '--log', dest='log_file', required=False, type=str,
                        default='cosmis_complex.log', help='''Logging file.''')
    parser.add_argument('--chain', dest='pdb_chain', default='A', type=str,
                        help='Chain ID of the subunit in the PBD file.')
    parser.add_argument('-d', dest='pdb_dir', required=True,
                        type=str, help='Directory containing the PDB files.')
    parser.add_argument('-w', '--overwrite', dest='overwrite', required=False,
                        action='store_true', help='''Whether to overwrite 
                        already computed MTR3D scores.''')
    parser.add_argument('-v', '--verbose', dest='verbose', required=False,
                        action='store_true', help='''Whether to output verbose
                        data: number of contacting residues and number of 
                        missense and synonymous variants in the neighborhood
                        of the mutation site.''')
    return parser.parse_args()



def get_ensembl_accession(record):
    """
    Convert a FASTA file record into an accession ID.

    Parameters
    ----------
    record : str
        Record ID is of the format: ">CCDS2.2|Hs109|chr1"

    Returns
    -------
    str
        Accession ID that can be used as dictionary key.

    """
    parts = record.id.split('.')
    return parts[0]


def get_uniprot_accession(record):
    """
    Convenience function to work with Biopython SeqIO.

    Parameters
    ----------
    record : Biopython SeqRecord
        A Biopython SeqRecord object.

    Returns
    -------
    str
        UniProt accession code.
    """
    parts = record.id.split('|')
    return parts[1]

def parse_config(config):
    """
    Parses the configuration file for computing COSMIS scores into a Python
    dictionary.

    Parameters
    ----------
    config : json
        Configuration file in JSON format.

    Returns
    -------
    dict
        A Python dictionary.

    """
    with open(config, 'rt') as ipf:
        configs = json.load(ipf)

    # do necessary sanity checks before return
    return configs


def count_variants(variants):
    """
    Collects the statistics about position-specific counts of missense and
    synonymous variants.

    Parameters
    ----------
    variants : list
        A list of variant identifiers: ['A123B', 'C456D']

    Returns
    -------
    dict
        A dictionary where the key is amino acid position and the value is
        the number of variants at this position. One dictionary for missense
        variants and one dictionary for synonymous variants.

    """
    missense_counts = defaultdict(int)
    synonymous_counts = defaultdict(int)
    for variant in variants:
        vv, ac, an = variant
        w = vv[0]  # wild-type amino acid
        v = vv[-1]  # mutant amino acid
        pos = vv[1:-1]  # position in the protein sequence
        # only consider rare variants
        if int(ac) / int(an) > 0.001:
            continue
        if w != v:  # missense variant
            missense_counts[int(pos)] += 1
        else:  # synonymous variant
            synonymous_counts[int(pos)] += 1
    return missense_counts, synonymous_counts


def get_transcript_info(
        uniprot_id, enst_ids, enst_ids_from_ensg, cds_dict, pep_dict, pep_dict2, variant_dict, enst_mp_counts
    ):
    """

    Parameters
    ----------
    transcript

    Returns
    -------

    """
    # Trying to find protein sequence of the uniprot ID
    try:
        pep_seq = pep_dict[uniprot_id]
    except KeyError:
        try:
            pep_seq = pep_dict2[uniprot_id]
        except KeyError:
            logging.critical(f'No peptide sequence found for {uniprot_id} in either pep_dict or pep_dict2.')
            # Handle the case where no peptide sequence is found
            #sys.exit(1)  # Optionally exit the program or handle the critical error appropriately

    # Empty list for storing enst ids
    valid_ensts = []

    # If no transcript ID was found from uniprot_to_enst
    if enst_ids is None:
        enst_ids = enst_ids_from_ensg       #assigning the new list to enst_ids
   
    for enst_id in enst_ids:
        try:
            cds_seq = cds_dict[enst_id].seq
        except KeyError:
            logging.critical(f'No CDS found for {enst_id} in Ensembl CDS database!')
            continue
        # skip if the CDS is incomplete
        if not seq_utils.is_valid_cds(cds_seq):
            print('Error: Invalid CDS.'.format(enst_id))
            continue
        if len(pep_seq) == len(cds_seq) // 3 - 1:
            valid_ensts.append(enst_id)
    if not valid_ensts:
        raise ValueError(
            'Error: {} are not compatible with {}.'.format(enst_ids, uniprot_id)
        )

    # get the one with most variable positions
    max_len = 0
    best_enst_id = valid_ensts[0]
    for enst_id in valid_ensts:
        try:
            var_pos = len(variant_dict[enst_id]['variants'])
        except KeyError:
            continue
        if max_len < var_pos:
            max_len = var_pos
            best_enst_id = enst_id

    # print(f'Best ENST ID for {uniprot_id} is {best_enst_id}.')

    # get the amino acid sequence of the transcript
    try:
        # Ensembl peptide ID for the transcript
        ensp_id = variant_dict[best_enst_id]['ensp'][0]
    except KeyError:
        logging.critical(f'Transcript {best_enst_id} not found in gnomAD database')
        #sys.exit(1)
        return None
        

    try:
        total_exp_mis_counts = enst_mp_counts[best_enst_id][-1]
        total_exp_syn_counts = enst_mp_counts[best_enst_id][-2]
    except KeyError:
        logging.critical(f'Transcript {best_enst_id} not found in transcript variant count file.')
        #sys.exit(1)
        return None
        

    # get all variants of this transcript reported in gnomAD
    try:
        variants = variant_dict[best_enst_id]['variants']
    except KeyError:
        logging.critical(f'No variants found for {best_enst_id} in gnomAD')
        #sys.exit(1)
        return None
        

    # get the coding sequence of the transcript
    try:
        transcript_cds = cds_dict[best_enst_id].seq
    except KeyError:
        logging.critical(f'No CDS found for {best_enst_id} in Ensembl CDS database!')
        transcript_cds = None
        return None

    # check that the CDS does not contain invalid nucleotides
    if not seq_utils.is_valid_cds(transcript_cds):
        print('ERROR: invalid CDS for', transcript_cds)
        #sys.exit(1)
        return None
        

    # get mutation rates and expected counts for each codon
    transcript_cds = transcript_cds[:-3]  # remove the stop codon
    codon_mutation_rates = seq_utils.get_codon_mutation_rates(transcript_cds)
    all_cds_ns_counts = seq_utils.count_poss_ns_variants(transcript_cds)

    # tabulate variants at each site
    # missense_counts and synonymous_counts are dictionary that maps
    # amino acid positions to variant counts
    missense_counts, synonymous_counts = count_variants(variants)

    # permutation test
    codon_mis_probs = [x[1] for x in codon_mutation_rates]
    codon_syn_probs = [x[0] for x in codon_mutation_rates]
    mis_p = codon_mis_probs / np.sum(codon_mis_probs)
    syn_p = codon_syn_probs / np.sum(codon_syn_probs)
    mis_pmt_matrix = seq_utils.permute_variants(
        total_exp_mis_counts, len(pep_seq), mis_p
    )
    syn_pmt_matrix = seq_utils.permute_variants(
        total_exp_syn_counts, len(pep_seq), syn_p
    )

    # return transcript statistics
    return ensp_id, codon_mutation_rates, missense_counts, \
           synonymous_counts, mis_pmt_matrix, syn_pmt_matrix, \
           all_cds_ns_counts, pep_seq


# New function to get ENST IDs from ENSG
def get_enst_from_ensg(ensg_id, configs):
    # Empty list for the ENST IDs
    enst_ids = []
    
    # Open the ENSEMBL CDS file
    with gzip.open(configs['ensembl_cds'], 'rt') as cds_handle:
        for record in SeqIO.parse(cds_handle, format='fasta'):
            # Extract the header
            header = record.description
            
            # Extract the ENSG ID from the header by splitting with ':'
            header_parts = header.split(':')
            if len(header_parts) > 6:
                extracted_ensg_id = header_parts[6].split('.')[0]  # Selecting the appropriate index and removing the period
                # Check if the extracted ENSG ID matches the specified ensg_id
                if extracted_ensg_id == ensg_id:
                    # Extract the ENST ID from the header by splitting with '>'
                    match = re.search(r'ENST(\d+)', header)
                    if match:
                        enst_id = 'ENST' + match.group(1)
                        enst_ids.append(enst_id)
    
    return enst_ids


#Preprocessing new fasta to make it work with SeqIO
def preprocess_fasta_header(line):
    # Function to remove spaces from the FASTA header line
    if line.startswith(">"):
        return ">" + line[1:].replace(" ", "")
    else:
        return line 

def load_datasets(configs):
    """

    Parameters
    ----------
    configs

    Returns
    -------
    """
    # ENSEMBL cds
    print('Reading ENSEMBL CDS database ...')
    with gzip.open(configs['ensembl_cds'], 'rt') as cds_handle:
        enst_cds_dict = SeqIO.to_dict(
            SeqIO.parse(cds_handle, format='fasta'),
            key_function=get_ensembl_accession
        )
    
    # ENSEMBL peptide sequences
    print('Reading UniProt protein sequence database ...')
    with gzip.open(configs['uniprot_pep'], 'rt') as pep_handle:
        pep_dict = SeqIO.to_dict(
            SeqIO.parse(pep_handle, format='fasta'),
            key_function=get_uniprot_accession
        )
    
    with gzip.open(configs['uniprot_pep2'], 'rt') as pep_handle:
        # Preprocess the file and create a string with modified headers
        processed_content = "\n".join(preprocess_fasta_header(line) for line in pep_handle)

    # Create a StringIO object to simulate a file-like object from the modified content
    processed_handle = StringIO(processed_content)

    # Now parse the modified content using SeqIO.parse()
    pep_dict2 = {}

    # Now parse the modified content using SeqIO.parse() and handle duplicate keys
    for record in SeqIO.parse(processed_handle, format='fasta'):
        try:
            key = get_uniprot_accession(record)
            if key not in pep_dict2:
                pep_dict2[key] = record
        except ValueError:
            # Ignore duplicates
            pass
    
        
    # parse gnomad transcript-level variants
    print('Reading gnomAD variant database ...')
    with open(configs['gnomad_variants'], 'rt') as variant_handle:
        # transcript_variants will be a dict of dicts where major version
        # ENSEMBL transcript IDs are the first level keys and "ccds", "ensp",
        # "swissprot", "variants" are the second level keys. The value of each
        # second-level key is a Python list.
        enst_variants = json.load(variant_handle)

    # Parse the file that maps Ensembl transcript IDs to PDB IDs
    with open(configs['uniprot_to_enst'], 'rt') as ipf:
        uniprot_to_enst = json.load(ipf)
    
    # Parse file that maps Ensembl gene ID to uniprot ID
    with open(configs['uniprot_to_ensg'], 'rt') as ipf:
        uniprot_to_ensg = json.load(ipf)

    # get transcript mutation probabilities and variant counts
    print('Reading transcript mutation probabilities and variant counts ...')
    enst_mp_counts = seq_utils.read_enst_mp_count(configs['enst_mp_counts'])

    return (enst_cds_dict, pep_dict, pep_dict2, enst_variants, uniprot_to_enst, uniprot_to_ensg, enst_mp_counts)


def main():
    
    # Parse command-line arguments
    args = parse_cmd()

    # Configure the logging system
    logging.basicConfig(
        filename=args.log_file,
        level=logging.INFO,
        filemode='w',
        format='%(levelname)s:%(asctime)s:%(message)s'
    )

    # Parse configuration file
    configs = parse_config(args.config)

    # Load datasets
    cds_dict, pep_dict, pep_dict2, variant_dict, uniprot_to_enst, uniprot_to_ensg, enst_mp_counts = load_datasets(configs)    
    
    # Storing all results from computation
    cosmis = []

    logging.info('Datasets have been loaded.')

    input_files = os.listdir(args.input)
    for input_file in input_files:
        with open(os.path.join(args.input, input_file), 'r') as ipf:
            # The first line of each file contains the pdb_file name
            pdb_file = ipf.readline().strip()
            chain_uniprot_pairs = {}
            for line in ipf.readlines():
                chain, uniprot_id = line.strip().split(':')
                chain_uniprot_pairs[chain] = uniprot_id

        logging.info(f'{input_file} being read')

        # uniprot ID of chain from argument  
        uniprot_id_for_chain = chain_uniprot_pairs[args.pdb_chain]

        # Empty dictionaries for storing values
        ensp_ids = {}
        codon_mutation_rates = {}
        mis_counts = {}
        syn_counts = {}
        mis_pmt_matrices = {}
        syn_pmt_matrices = {}
        all_cds_ns_counts = {}
        pep_seqs = {}


        # Flag to break the loop later if needed
        found = False

        # Looping through all chain: Uniprot ID pairs in the txt files 
        for c, uniprot_id in chain_uniprot_pairs.items():
            logging.info(f'Processing {c},{uniprot_id} from {input_file}')
    
            # Trying to find transcript ID (ENST)
            try:
                enst_ids = uniprot_to_enst[uniprot_id] #checking first dataset. Gets ENST IDs from Uniprot ID
                ensg_id = None
            
            except KeyError:            #if KeyError detected
                logging.critical(
                        f'No transcript IDs were mapped to {uniprot_id} in uniprot_to_enst'
                    )
                
                # Setting enst id as None if its not found
                enst_ids = None 
                
                # Trying to find gene ID instead
                try:
                    ensg_id = uniprot_to_ensg[uniprot_id]
                except KeyError:
                    logging.critical(f'No Ensembl gene IDs were mapped to {uniprot_id} in uniprot_to_ensg.')

                    # Activate flag and break
                    found = True
                    break
            
            # Extract the ENST IDs from ENSG IDs
            if enst_ids == None:
                enst_ids_from_ensg = get_enst_from_ensg(ensg_id, configs)
            else:
                enst_ids_from_ensg = None


            try:
                results = get_transcript_info(
            uniprot_id, enst_ids, enst_ids_from_ensg, cds_dict, pep_dict, pep_dict2, variant_dict, enst_mp_counts
        )
                if results == None:
                    found=True
                    break
            except ValueError:
                logging.critical('No valid CDS found for {}.'.format(uniprot_id))
                #sys.exit(1)
                found = True
                break
            except KeyError:
                logging.critical('No transcript record found for {} in gnomAD.'.format(uniprot_id))
                #sys.exit(1)
                found = True
                break

              
               
            logging.info(f'ENST processed for: {uniprot_id}, {results[0]}, {pdb_file}')

            ensp_ids[c] = results[0]
            codon_mutation_rates[c] = results[1]
            mis_counts[c] = results[2]
            syn_counts[c] = results[3]
            mis_pmt_matrices[c] = results[4]
            syn_pmt_matrices[c] = results[5]
            all_cds_ns_counts[c] = results[6]
            pep_seqs[c] = results[7]


            logging.info(f'{c},{uniprot_id} stored in dictionaries!')

            # index all contacts by residue ID
        if found:
            found=False
            continue  
            
        pdb_dir = args.pdb_dir  # Get the directory containing the PDB files
        pdb_file_path = os.path.join(pdb_dir, pdb_file) 
        pdb_parser = PDBParser(PERMISSIVE=1)
        structure = pdb_parser.get_structure(id='NA', file=pdb_file_path)
        all_aa_residues = [aa for aa in structure[0].get_residues() if is_aa(aa)]
        all_contacts = pdb_utils.search_for_all_contacts(all_aa_residues, radius=8)
        indexed_contacts = defaultdict(list)

        logging.info(f'PDB parsed for {input_file}')

        for c in all_contacts:
            indexed_contacts[c.get_res_a()].append(c.get_res_b())
            indexed_contacts[c.get_res_b()].append(c.get_res_a())

        chain_id = args.pdb_chain
        
        chain_struct = structure[0][chain_id]
        for seq_pos, seq_aa in enumerate(pep_seqs[chain_id], start=1):
            # check that the amino acid in ENSP sequence matches
            # that in the PDB structure
            try:
                res = chain_struct[seq_pos]
            except KeyError:
                print(
                    'Residue %s not found in chain %s in PDB file: %s' %
                    (seq_pos, chain_id, pdb_file_path)
                )
                logging.info(f'Residue {seq_pos} not found in {chain_id} in pdb file {pdb_file_path}')
                continue
            pdb_aa = seq1(res.get_resname())
            if seq_aa != pdb_aa:
                print('Residue in ENSP did not match that in PDB at', seq_pos)
                logging.info(f'Residue in ENSP did not match that in PDB at {seq_pos}')
                continue
            contact_res = indexed_contacts[res]
            print('Current residue:', res.get_full_id())
            total_missense_obs = mis_counts[chain_id].setdefault(seq_pos, 0)
            total_synonymous_obs = syn_counts[chain_id].setdefault(seq_pos, 0)
            total_missense_poss = all_cds_ns_counts[chain_id][seq_pos - 1][0]
            total_synonyms_poss = all_cds_ns_counts[chain_id][seq_pos - 1][1]
            total_synonymous_rate = codon_mutation_rates[chain_id][seq_pos - 1][0]
            total_missense_rate = codon_mutation_rates[chain_id][seq_pos - 1][1]
            flag=False
            for c_res in contact_res:
                # count the total # observed variants in contacting residues
                print('\tContacts:', c_res.get_full_id())
                contact_chain_id = c_res.get_full_id()[2]
                if contact_chain_id not in mis_counts or contact_chain_id not in syn_counts:
                    logging.info(f"Chain {contact_chain_id} not found in mis_counts or syn_counts. breaking...")
                    flag=True
                    break
                j = c_res.get_full_id()[3][1]
                total_missense_obs += mis_counts[contact_chain_id].setdefault(j, 0)
                total_synonymous_obs += syn_counts[contact_chain_id].setdefault(j, 0)

                # count the total # expected variants
                try:
                    total_missense_poss += all_cds_ns_counts[contact_chain_id][j - 1][0]
                    total_synonyms_poss += all_cds_ns_counts[contact_chain_id][j - 1][1]
                    total_synonymous_rate += codon_mutation_rates[contact_chain_id][j - 1][0]
                    total_missense_rate += codon_mutation_rates[contact_chain_id][j - 1][1]
                except IndexError:
                    logging.info('{} not in CDS that has {} residues.'.format(j, len(
                        all_cds_ns_counts)))
                    #sys.exit(1)
                    flag=True
                    break
            if flag:
                flag=False
                break
            mis_pmt_mean, mis_pmt_sd, mis_p_value = seq_utils.get_permutation_stats(
                mis_pmt_matrices, contact_res + [res], total_missense_obs
            )
            syn_pmt_mean, syn_pmt_sd, syn_p_value = seq_utils.get_permutation_stats(
                syn_pmt_matrices, contact_res + [res], total_synonymous_obs
            )


            # add entries to be written to disk file
            cosmis.append(
                [
                    uniprot_id_for_chain, seq_pos, seq_aa,
                    total_synonyms_poss,
                    total_missense_poss,
                    '%.3e' % total_synonymous_rate,
                    total_synonymous_obs,
                    '%.3e' % total_missense_rate,
                    total_missense_obs,
                    '{:.3f}'.format(mis_pmt_mean),
                    '{:.3f}'.format(mis_pmt_sd),
                    '{:.3e}'.format(mis_p_value),
                    '{:.3f}'.format(syn_pmt_mean),
                    '{:.3f}'.format(syn_pmt_sd),
                    '{:.3e}'.format(syn_p_value),
                ]
            )

        logging.info(f'{input_file} info added')

    with open(file=args.output_file, mode='wt') as opf:
        header = [
                        'uniprot_id', 'ensp_pos', 'ensp_aa',
                        'cs_syn_poss', 'cs_mis_poss',
                        'cs_syn_prob', 'cs_syn_obs', 'cs_mis_prob', 'cs_mis_obs',
                        'mis_pmt_mean', 'mis_pmt_sd', 'mis_p_value', 'syn_pmt_mean',
                        'syn_pmt_sd', 'syn_p_value'
                    ]
        csv_writer = csv.writer(opf, delimiter='\t')
        csv_writer.writerow(header)
        csv_writer.writerows(cosmis)
        
    logging.info('output file created!')

if __name__ == '__main__':
    main()
