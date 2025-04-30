#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=23:00:00
#SBATCH --job-name=z1_seperatingFireboxGalaxies
#SBATCH --output=z1_seperatingFireboxGalaxies.out
#SBATCH --error=z1_seperatingFireboxGalaxies.err

# module purge 
# ml python/3.11.5

# cd /scratch/m/murray/dtolgay/post_processing_fire_outputs/firebox_halo_finder

# # Run the virtual environment
# source venv/bin/activate

# export PYTHONUNBUFFERED=1 # this is the write .out file while the code is running 
# python seperate_firebox_galaxies_from_largeBox.py 2.0 900 1000



module purge 
ml python/3.11.5

cd /scratch/m/murray/dtolgay/post_processing_fire_outputs/firebox_halo_finder

source venv/bin/activate

export PYTHONUNBUFFERED=1
python seperate_firebox_galaxies_from_largeBox.py 1.0 900 1000