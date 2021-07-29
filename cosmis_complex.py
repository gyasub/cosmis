#!/usr/bin/env python3

"""
    @summary: This is the single-protein version of the COSMIS method, hence its
    name 'cosmis_sp.py'. It was designed to work with proteins for which the
    COSMIS scores were not already computed, because the protein structure
    database used wasn't up to date. However, if there exists a structure or
    homology model for the protein of interest, its COSMIS scores can still be
    computed using this script.
    This script assumes that residues in the PDB file are numbered according to
    the amino acid sequence in UniProt.
    @author: Bian Li
    @contact: bian.li@vanderbilt.edu
    @change: Last modified 2/12/2021.

"""

import csv
import gzip
import json
import sys
import numpy as np
from argparse import ArgumentParser
from collections import defaultdict

from Bio import SeqIO
from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils import seq1

from cosmis.utils import pdb_utils, seq_utils


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
                        type=str, help='''Pairs of PDB chain ID:Ensembl 
                        transcript ID associated with the input PDB file.''')
    parser.add_argument('-p', '--pdb', dest='pdb_file', required=True,
                        type=str, help='''PDB file containing the structure
                        of the protein structure of the given transcript.''')
    parser.add_argument('-o', '--output', dest='output_file', required=True,
                        type=str, help='''Output file to store the COSMIS scores
                        of the protein.''')
    parser.add_argument('--chain', dest='pdb_chain', default='A', type=str,
                        help='Chain ID of the subunit in the PBD file.')
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


def get_ccds_accession(record):
    """

    Parameters
    ----------
    record

    Returns
    -------

    """
    parts = record.id.split('|')
    return parts[0]


def get_transcript_pep_seq(enst_id, ensp_id, pep_dict):
    """

    Parameters
    ----------
    enst_id : str

    ensp_id : str

    pep_dict : dict

    Returns
    -------

    """
    try:
        transcript_pep = pep_dict[ensp_id].seq
    except KeyError:
        print('%s not found in given database' % ensp_id)
        print('%s was skipped ...' % enst_id)
        return None
    return transcript_pep


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
    #
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


def parse_ensembl_cds(file=None):
    """

    Parameters
    ----------
    file

    Returns
    -------

    """
    print('Parsing ENSEMBL concensus coding sequence database ...')
    with gzip.open(file, mode='rt') as cds_handle:
        ensembl_cds_dict = SeqIO.to_dict(
            SeqIO.parse(cds_handle, format='fasta'),
            key_function=get_ensembl_accession
        )
    return ensembl_cds_dict


def parse_ccds(file=None):
    """

    Parameters
    ----------
    file

    Returns
    -------

    """
    print('Parsing NCBI CCDS database ...')
    with gzip.open(file, mode='rt') as ccds_handle:
        ccds_dict = SeqIO.to_dict(
            SeqIO.parse(ccds_handle, format='fasta'),
            key_function=get_ccds_accession
        )
    return ccds_dict


def parse_ensembl_pep(file=None):
    """

    Parameters
    ----------
    file

    Returns
    -------

    """
    print('Parsing Ensembl protein sequence database ...')
    with gzip.open(file, mode='rt') as pep_handle:
        pep_dict = SeqIO.to_dict(
            SeqIO.parse(pep_handle, format='fasta'),
            key_function=get_ensembl_accession
        )
    return pep_dict


def get_transcript_info(
        transcript, ensembl_cds_dict, ccds_dict, pep_dict, gnomad_variants, enst_mp_counts
    ):
    """

    Parameters
    ----------
    transcript

    Returns
    -------

    """
    # get the amino acid sequence of the transcript
    try:
        # Ensembl peptide ID for the transcript
        ensp_id = gnomad_variants[transcript]['ensp'][0]
    except KeyError:
        print(
            'Transcript %s not found in gnomAD database',
            transcript
        )
        sys.exit(1)

    try:
        total_exp_mis_counts = enst_mp_counts[transcript][-1]
        total_exp_syn_counts = enst_mp_counts[transcript][-2]
    except KeyError:
        print(
            'Transcript {} not found.'.format(transcript)
        )
        sys.exit(1)

    # get the peptide sequence from peptide sequence database
    transcript_pep_seq = get_transcript_pep_seq(
        transcript, ensp_id, pep_dict
    )

    # get all variants of this transcript reported in gnomAD
    try:
        variants = gnomad_variants[transcript]['variants']
    except KeyError:
        print('No variants found for %s in gnomAD', transcript)
        sys.exit(1)

    # get the coding sequence of the transcript
    try:
        transcript_cds = ensembl_cds_dict[transcript].seq
    except KeyError:
        print(
            '''No CDS found in Ensembl CDS database! 
            Looking for it in the CCDS database ...'''
        )
        transcript_cds = None

    if transcript_cds is None:
        try:
            ccds_id = gnomad_variants[transcript]['ccds'][0]
            transcript_cds = ccds_dict[ccds_id].seq
        except KeyError:
            print('ERROR: No CDS found in CCDS database!')
            sys.exit(1)

    # check that the CDS does not contain invalid nucleotides
    if not seq_utils.is_valid_cds(transcript_cds):
        print('ERROR: invalid CDS for', transcript_cds)
        sys.exit(1)

    # get peptide sequence by translating the CDS
    # if peptide sequence isn't in Ensembl
    if transcript_pep_seq is None:
        transcript_pep_seq = seq_utils.translate(transcript_cds)

    uniprot_ids = gnomad_variants[transcript]['swissprot']
    if len(uniprot_ids) > 1:
        print(
            'ERROR: more than one UniProt IDs found for transcript '
            '%s: %s', transcript, ','.join(uniprot_ids)
        )
        sys.exit(1)

    # get the UniProt ID
    uniprot_id = uniprot_ids[0]
    if not uniprot_id:
        uniprot_id = "NA"  # UniProt IDs may be missing in gnomAD

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
        total_exp_mis_counts, len(transcript_pep_seq), mis_p
    )
    syn_pmt_matrix = seq_utils.permute_variants(
        total_exp_syn_counts, len(transcript_pep_seq), syn_p
    )

    # return transcript statistics
    return ensp_id, uniprot_id, codon_mutation_rates, missense_counts, \
           synonymous_counts, mis_pmt_matrix, syn_pmt_matrix, \
           all_cds_ns_counts, transcript_pep_seq


def main():
    # parse command-line arguments
    args = parse_cmd()

    # parse configuration file
    configs = parse_config(args.config)

    # ENSEMBL cds
    ensembl_cds_dict = parse_ensembl_cds(configs['ensembl_cds'])
    ccds_dict = parse_ccds(configs['ccds_cds'])
    pep_dict = parse_ensembl_pep(configs['ensembl_pep'])

    # parse gnomad transcript-level variants
    print('Reading gnomAD variant database ...')
    with open(configs['gnomad_variants'], 'rt') as variant_handle:
        # transcript_variants will be a dict of dicts where major version
        # ENSEMBL transcript IDs are the first level keys and "ccds", "ensp",
        # "swissprot", "variants" are the second level keys. The value of each
        # second-level key is a Python list.
        transcript_variants = json.load(variant_handle)

    # get transcript mutation probabilities and variant counts
    print('Reading transcript mutation probabilities and variant counts ...')
    enst_mp_counts = seq_utils.read_enst_mp_count(configs['enst_mp_counts'])

    # parse chain-transcript pairs
    with open(args.input, 'rt') as ipf:
        chain_transcript_pairs = {
            c: t for c, t in [line.strip().split(':') for line in ipf]
        }
    transcript = chain_transcript_pairs[args.pdb_chain]

    #
    cosmis_scores = []
    ensp_ids = {}
    uniprot_ids = {}
    codon_mutation_rates = {}
    mis_counts = {}
    syn_counts = {}
    mis_pmt_matrices = {}
    syn_pmt_matrices = {}
    all_cds_ns_counts = {}
    transcript_pep_seqs = {}
    for c, t in chain_transcript_pairs.items():
        results = get_transcript_info(
            t, ensembl_cds_dict, ccds_dict, pep_dict, transcript_variants, enst_mp_counts
        )
        ensp_ids[c] = results[0]
        uniprot_ids[c] = results[1]
        codon_mutation_rates[c] = results[2]
        mis_counts[c] = results[3]
        syn_counts[c] = results[4]
        mis_pmt_matrices[c] = results[5]
        syn_pmt_matrices[c] = results[6]
        all_cds_ns_counts[c] = results[7]
        transcript_pep_seqs[c] = results[8]
    # index all contacts by residue ID
    pdb_parser = PDBParser(PERMISSIVE=1)
    structure = pdb_parser.get_structure(id='NA', file=args.pdb_file)
    all_aa_residues = [aa for aa in structure[0].get_residues() if is_aa(aa)]
    all_contacts = pdb_utils.search_for_all_contacts(all_aa_residues, radius=8)
    indexed_contacts = defaultdict(list)
    for c in all_contacts:
        indexed_contacts[c.get_res_a()].append(c.get_res_b())
        indexed_contacts[c.get_res_b()].append(c.get_res_a())

    chain_id = args.pdb_chain
    chain_struct = structure[0][chain_id]
    for seq_pos, seq_aa in enumerate(transcript_pep_seqs[chain_id], start=1):
        # check that the amino acid in ENSP sequence matches
        # that in the PDB structure
        try:
            res = chain_struct[seq_pos]
        except KeyError:
            print(
                'Residue %s not found in chain %s in PDB file: %s' %
                (seq_pos, chain_id, args.pdb_file)
            )
            continue
        pdb_aa = seq1(res.get_resname())
        if seq_aa != pdb_aa:
            print('Residue in ENSP did not match that in PDB at', seq_pos)
            continue
        contact_res = indexed_contacts[res]
        print('Current residue:', res.get_full_id())
        total_missense_obs = mis_counts[chain_id].setdefault(seq_pos, 0)
        total_synonymous_obs = syn_counts[chain_id].setdefault(seq_pos, 0)
        total_missense_poss = all_cds_ns_counts[chain_id][seq_pos - 1][0]
        total_synonyms_poss = all_cds_ns_counts[chain_id][seq_pos - 1][1]
        total_synonymous_rate = codon_mutation_rates[chain_id][seq_pos - 1][0]
        total_missense_rate = codon_mutation_rates[chain_id][seq_pos - 1][1]
        for c_res in contact_res:
            # count the total # observed variants in contacting residues
            print('\tContacts:', c_res.get_full_id())
            contact_chain_id = c_res.get_full_id()[2]
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
                print('{} not in CDS that has {} residues.'.format(j, len(
                    all_cds_ns_counts)))
                sys.exit(1)
        mis_pmt_mean, mis_pmt_sd, mis_p_value = seq_utils.get_permutation_stats(
            mis_pmt_matrices, contact_res + [res], total_missense_obs
        )
        syn_pmt_mean, syn_pmt_sd, syn_p_value = seq_utils.get_permutation_stats(
            syn_pmt_matrices, contact_res + [res], total_synonymous_obs
        )

        # add entries to be written to disk file
        cosmis_scores.append(
            [
                transcript, seq_pos, seq_aa,
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

    with open(file=args.output_file, mode='wt') as opf:
        header = [
            'enst_id', 'ensp_pos', 'ensp_aa',
            'cs_syn_poss', 'cs_mis_poss',
            'cs_syn_prob', 'cs_syn_obs', 'cs_mis_prob', 'cs_mis_obs',
            'mis_pmt_mean', 'mis_pmt_sd', 'mis_p_value', 'syn_pmt_mean',
            'syn_pmt_sd', 'syn_p_value'
        ]
        csv_writer = csv.writer(opf, delimiter='\t')
        csv_writer.writerow(header)
        csv_writer.writerows(cosmis_scores)


if __name__ == '__main__':
    main()