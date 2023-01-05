#!/bin/bash

#SBATCH --job-name=BASECASE_5yr_CY
#SBATCH --time=0-25:00:00
#SBATCH -A sera
#SBATCH --mem=150G
#SBATCH -o BASECASE_5yr_CY.log
#SBATCH --mail-user=mark.chung@nrel.gov
#SBATCH --mail-type=ALL
module purge
module load conda
#export PATH=/nopt/nrel/apps/xpressmp/8.8.0/lib:$PATH
module load xpressmp/8.8.0
module load gurobi/9.0.2

#alias julia='/home/mchung/julia_software/julia-1.8.4/bin/julia'
export PATH=$PATH:/lustre/eaglefs/projects/mbap/sera/mchung/julia-1.8.4/bin/

cd /lustre/eaglefs/projects/mbap/sera/vvijayak/SERA.jl

julia --project /lustre/eaglefs/projects/mbap/sera/mchung/supply-chain-inputs-external/SERA-supply-chain-inputs/cases/BASECASE_5yr_CY/run_case.jl
# julia --project /lustre/eaglefs/projects/mbap/sera/vvijayak/supply-chain-inputs/cases/BASECASE_Centralonly-15util-maxcap-rolling/run_case.jl

# export LD_LIBRARY_PATH=/projects/mbap/sera/bin/lib



