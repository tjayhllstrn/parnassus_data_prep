#!/bin/bash

#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem-per-cpu=5000
#SBATCH --nodes=2
#SBATCH --chdir=/scratch/slurm
#SBATCH --output=/farm_out/%u/%x-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%j-%N.err
#SBATCH --job-name=pipeline0_runSingleHipo_pippi0_RGAMCpass2
#SBATCH --time=5:00:00

cd /w/hallb-scshelf2102/clas12/users/tjhellst/parnassus_data_prep

root -l -b -q 'macros/hipo2tree_parnassus.C("/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb/clasdis_rga_fa18_inb_45nA_10604MeV-0001.hipo","prepped_data/clasdis_rga_fa18_inb_45nA_10604MeV-0001.root",-1)'
