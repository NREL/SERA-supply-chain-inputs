#!/bin/bash

#SBATCH --job-name=BASECASE_Centralonly-lowercap_station_plant-pipeline_e5
#SBATCH --time=0-30:00:00
#SBATCH -A sera
#SBATCH --mem=150G
#SBATCH -o BASECASE_Centralonly-lowercap_station_plant-pipeline_e5.log
#SBATCH --mail-user=vishnu.vijayakumar@nrel.gov
#SBATCH --mail-type=ALL
module purge
#export PATH=/nopt/nrel/apps/xpressmp/8.8.0/lib:$PATH
module load xpressmp/8.8.0
module load conda
conda activate /projects/mbap/sera/vvijayak/SERA.jl/julia_software
alias julia='julia'

cd /lustre/eaglefs/projects/mbap/sera/vvijayak/supply-chain-inputs

julia --project /lustre/eaglefs/projects/mbap/sera/vvijayak/supply-chain-inputs/cases/BASECASE_Centralonly-lowercap_station_plant-pipeline_e5/run_case.jl

export LD_LIBRARY_PATH=/projects/mbap/sera/bin/lib



