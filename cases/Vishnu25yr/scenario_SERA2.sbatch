#!/bin/bash
#SBATCH --mem=150G
#SBATCH --job-name=Vishnu25yr
#SBATCH --time=0-25:00:00
#SBATCH -A sera
#SBATCH --chdir=/home/ccyang/SERA-supply-chain-inputs/cases/Vishnu25yr
#SBATCH --output=Vishnu25yr.out 
#SBATCH --error=Vishnu25yr.err
#SBATCH --mail-user=ccyang@ucdavis.edu
#SBATCH --mail-type=ALL
module purge
#export PATH=/nopt/nrel/apps/xpressmp/8.8.0/lib:$PATH
module load xpressmp/8.8.0
module load gurobi/9.0.2
module load conda

cd /home/ccyang/SERA-supply-chain-inputs/cases/Vishnu25yr/

julia --project /home/ccyang/SERA-supply-chain-inputs/cases/Vishnu25yr/run_case.jl

# export LD_LIBRARY_PATH=/projects/mbap/sera/bin/lib



