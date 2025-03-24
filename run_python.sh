#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=40
#SBATCH --time=23:00:00
#SBATCH --job-name=5nodes_z0_seperatingFireboxGalaxies
#SBATCH --output=5nodes_z0_seperatingFireboxGalaxies.out
#SBATCH --error=5nodes_z0_seperatingFireboxGalaxies.err

module purge 
ml python/3.11.5

cd /scratch/m/murray/dtolgay/post_processing_fire_outputs/firebox_halo_finder

# Run the virtual environment
source venv/bin/activate

export PYTHONUNBUFFERED=1 # this is the write .out file while the code is running 
python seperate_firebox_galaxies_from_largeBox.py 0.0