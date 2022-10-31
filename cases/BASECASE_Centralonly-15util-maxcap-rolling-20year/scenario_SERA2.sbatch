#!/bin/bash
#SBATCH --mem=150G
#SBATCH --job-name=BASECASE_Centralonly-15util-maxcap-rolling-20year
#SBATCH --time=0-25:00:00
#SBATCH -A mbap
#SBATCH -o BASECASE_Centralonly-15util-maxcap-rolling-25year_30tpd_maxgh2trucklink3.log
#SBATCH --mail-user=vishnu.vijayakumar@nrel.gov
#SBATCH --mail-type=ALL
module purge
#export PATH=/nopt/nrel/apps/xpressmp/8.8.0/lib:$PATH
module load xpressmp/8.8.0
module load gurobi/9.0.2
module load conda
conda activate /projects/mbap/sera/vvijayak/SERA.jl/julia_software
alias julia='julia'

cd /lustre/eaglefs/projects/mbap/sera/vvijayak/carbon-neutral-california

julia --project /lustre/eaglefs/projects/mbap/sera/vvijayak/carbon-neutral-california/cases/BASECASE_Centralonly-15util-maxcap-rolling-20year/run_case.jl

export LD_LIBRARY_PATH=/projects/mbap/sera/bin/lib



