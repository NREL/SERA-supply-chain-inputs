#!/bin/bash

#SBATCH --job-name=BASECASE_Centralonly-15util-maxcap-rolling-50percent
#SBATCH --time=0-25:00:00
#SBATCH -A mbap
#SBATCH --mem=150G
#SBATCH -o BASECASE_Centralonly-15util-maxcap-rolling-50percent_latest_gh2truckterminal_link_30tpdmaxcap.log
#SBATCH --mail-user=vishnu.vijayakumar@nrel.gov
#SBATCH --mail-type=ALL
module purge
#export PATH=/nopt/nrel/apps/xpressmp/8.8.0/lib:$PATH
module load xpressmp/8.8.0
module load conda
conda activate /projects/mbap/sera/vvijayak/SERA.jl/julia_software
alias julia='julia'

cd /lustre/eaglefs/projects/mbap/sera/vvijayak/supply-chain-inputs

julia --project /lustre/eaglefs/projects/mbap/sera/vvijayak/supply-chain-inputs/cases/BASECASE_Centralonly-15util-maxcap-rolling-50percent/run_case.jl

export LD_LIBRARY_PATH=/projects/mbap/sera/bin/lib



