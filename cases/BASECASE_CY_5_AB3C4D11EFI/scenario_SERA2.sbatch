#!/bin/bash

#SBATCH --job-name=CY_5_AB3C4D11EFI
#SBATCH --time=0-30:00:00
#SBATCH -A sera
#SBATCH --output=BASECASE_CY_5_AB3C4D11EFI.out 
#SBATCH --error=BASECASE_CY_5_AB3C4D11EFI.err
#SBATCH --mail-user=ccyang@ucdavis.edu
#SBATCH --mail-type=ALL
module purge
module load conda
export PATH=/nopt/nrel/apps/xpressmp/8.8.0/lib:$PATH
module load xpressmp/8.8.0

#alias julia='/home/ccyang/julia_software/julia-1.7.3/bin/julia'
export PATH=$PATH:/home/ccyang/julia_software/julia-1.7.3/bin/

cd /home/ccyang/SERA-supply-chain-inputs/cases/BASECASE_CY_5_AB3C4D11EFI/

julia --project /home/ccyang/SERA-supply-chain-inputs/cases/BASECASE_CY_5_AB3C4D11EFI/run_case.jl
# julia --project /lustre/eaglefs/projects/mbap/sera/vvijayak/supply-chain-inputs/cases/BASECASE_Centralonly-15util-maxcap-rolling/run_case.jl

# export LD_LIBRARY_PATH=/projects/mbap/sera/bin/lib



