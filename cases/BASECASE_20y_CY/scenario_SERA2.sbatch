#!/bin/bash
#SBATCH --job-name=BASECASE_20y_CY
#SBATCH --time=0-25:00:00
#SBATCH -A sera
#SBATCH -o BASECASE_20y_CY.log
#SBATCH --mail-user=ccyang@ucdavis.edu
#SBATCH --mail-type=ALL
module purge
#export PATH=/nopt/nrel/apps/xpressmp/8.8.0/lib:$PATH
module load xpressmp/8.8.0
module load gurobi/9.0.2
module load conda
conda activate /projects/mbap/sera/vvijayak/SERA.jl/julia_software
export PATH=$PATH:/home/ccyang/julia_software/julia-1.7.3/bin/

cd /lustre/eaglefs/projects/sera/SERA2.0/carbon-neutral-california/SERA.jl

julia --project /home/ccyang/SERA-supply-chain-inputs/cases/BASECASE_20yr_CY/run_case.jl

# export LD_LIBRARY_PATH=/projects/mbap/sera/bin/lib



