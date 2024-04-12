#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l mem_free=32G
#$ -l gpu_mem=24000
#$ -l scratch=24G
#$ -l h_rt=30:00:00
#$ -l compute_cap=61
#$ -q gpu.q


date
hostname

module load CBI miniconda3/23.5.2-0-py311
conda activate cosmis

python scripts/cosmis_complex_batch_final.py -c data_paths.json -i /wynton/home/rotation/gyasu/FoldDock/output_txt_dir/ -o AFcomplex_chainA_run2.tsv -l AFcomplex_chainA_run2.log --chain A -d /wynton/home/rotation/gyasu/FoldDock/output_pdbs/


### qsub /wynton/home/rotation/gyasu/capra_lab/cosmis/run_cosmis/run_cosmis.sh

#historical commands
#python scripts/cosmis_complex_batch_final.py -c data_paths.json -i /wynton/home/rotation/gyasu/FoldDock/debug_txt_dir/ -o debug_test.tsv -l debug.log --chain A -d /wynton/home/rotation/gyasu/FoldDock/output_pdbs/
#python scripts/cosmis_complex_batch_final.py -c data_paths.json -i /wynton/home/rotation/gyasu/FoldDock/output_txt_dir/ -o AFcomplex_chainB.tsv -l AFcomplex_chainB.log --chain B -d /wynton/home/rotation/gyasu/FoldDock/output_pdbs/
