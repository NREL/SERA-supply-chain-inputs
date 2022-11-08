#!/bin/bash
#SBATCH --mem=150G
#SBATCH --job-name=NS1_Vishnu25yr_3
#SBATCH --time=0-25:00:00
#SBATCH -A sera
#SBATCH --chdir=/home/ccyang/SERA-supply-chain-inputs/cases/NS1_Vishnu25yr_3
#SBATCH --output=NS1_Vishnu25yr_3.out 
#SBATCH --error=NS1_Vishnu25yr_3.err
#SBATCH --mail-user=ccyang@ucdavis.edu
#SBATCH --mail-type=ALL
module purge
#export PATH=/nopt/nrel/apps/xpressmp/8.8.0/lib:$PATH
module load xpressmp/8.8.0
module load gurobi/9.0.2
module load conda

cd /home/ccyang/SERA-supply-chain-inputs/cases/NS1_Vishnu25yr_3/

julia --project /home/ccyang/SERA-supply-chain-inputs/cases/NS1_Vishnu25yr_3/run_case.jl

# export LD_LIBRARY_PATH=/projects/mbap/sera/bin/lib



