#!/bin/bash

#SBATCH --job-name=CY_15_AB3C4D8EF
#SBATCH --time=0-65:00:00
#SBATCH -A sera
#SBATCH --mem=150G
#SBATCH -o BASECASE_CY_15_AB3C4D8EF.log
#SBATCH --mail-user=ccyang@ucdavis.edu
#SBATCH --mail-type=ALL
module purge
#export PATH=/nopt/nrel/apps/xpressmp/8.8.0/lib:$PATH
module load xpressmp/8.8.0
module load gurobi/9.0.2
module load conda

#alias julia='/home/ccyang/julia_software/julia-1.7.3/bin/julia'
export PATH=$PATH:/home/ccyang/julia_software/julia-1.7.3/bin/

cd /lustre/eaglefs/projects/sera/SERA2.0/carbon-neutral-california/SERA.jl

julia --project /home/ccyang/SERA-supply-chain-inputs/cases/BASECASE_CY_15_AB3C4D8EF/run_case.jl
# julia --project /lustre/eaglefs/projects/mbap/sera/vvijayak/supply-chain-inputs/cases/BASECASE_Centralonly-15util-maxcap-rolling/run_case.jl

# export LD_LIBRARY_PATH=/projects/mbap/sera/bin/lib



