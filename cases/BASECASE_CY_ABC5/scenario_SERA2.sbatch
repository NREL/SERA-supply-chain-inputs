#!/bin/bash

#SBATCH --job-name=BASECASE_CY_AB5
#SBATCH --time=0-25:00:00
#SBATCH -A sera
#SBATCH -o BASECASE_CY_AB5.log
#SBATCH --mail-user=ccyang@ucdavis.edu
#SBATCH --mail-type=ALL
module purge
module load conda
export PATH=/nopt/nrel/apps/xpressmp/8.8.0/lib:$PATH
module load xpressmp/8.8.0

#alias julia='/home/ccyang/julia_software/julia-1.7.3/bin/julia'
export PATH=$PATH:/home/ccyang/julia_software/julia-1.7.3/bin/

cd /lustre/eaglefs/projects/sera/SERA2.0/carbon-neutral-california/SERA.jl

julia --project /home/ccyang/SERA-supply-chain-inputs/cases/BASECASE_CY_AB5/run_case.jl
# julia --project /lustre/eaglefs/projects/mbap/sera/vvijayak/supply-chain-inputs/cases/BASECASE_Centralonly-15util-maxcap-rolling/run_case.jl

# export LD_LIBRARY_PATH=/projects/mbap/sera/bin/lib



