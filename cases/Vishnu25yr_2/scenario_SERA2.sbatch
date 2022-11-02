#!/bin/bash
#SBATCH --mem=150G
#SBATCH --job-name=Vishnu25yr
#SBATCH --time=0-55:00:00
#SBATCH -A sera
#SBATCH --chdir=/home/ccyang/SERA-supply-chain-inputs/cases/BASECASE_CY_5_AB3C4D10EFG3
#SBATCH -o Vishnu25yr_2.out
#SBATCH --error=Vishnu25yr_2.err
#SBATCH --mail-user=ccyang@ucdavis.edu
#SBATCH --mail-type=ALL
module purge
#export PATH=/nopt/nrel/apps/xpressmp/8.8.0/lib:$PATH
module load xpressmp/8.8.0
module load gurobi/9.0.2
module load conda

cd /home/ccyang/SERA-supply-chain-inputs/cases/Vishnu25yr_2/

julia --project /home/ccyang/SERA-supply-chain-inputs/cases/Vishnu25yr_2/run_case.jl

# export LD_LIBRARY_PATH=/projects/mbap/sera/bin/lib



