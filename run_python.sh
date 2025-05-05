#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=23:00:00
#SBATCH --job-name=__NAME__
#SBATCH --output=__NAME__.out
#SBATCH --error=__NAME__.err

module purge 
ml python/3.11.5

cd /scratch/m/murray/dtolgay/post_processing_fire_outputs/firebox_halo_finder
source venv/bin/activate
export PYTHONUNBUFFERED=1

# Parameters
z=__Z__
start=__START__
end=__END__
n_galaxies=10

# Calculate chunk size
step=$(( (end - start) / n_galaxies ))

# Run in background
for ((i=0; i<n_galaxies; i++)); do
    gal_start=$((start + i * step))

    # Use the last galaxy number as the end of the final chunk
    if [ $i -eq $((n_galaxies - 1)) ]; then
        gal_end=$end
    else
        gal_end=$((gal_start + step))
    fi

    echo "Processing galaxies $gal_start to $gal_end"
    python seperate_firebox_galaxies_from_largeBox.py $z $gal_start $gal_end &
done

# Wait for all background jobs to finish
wait
echo "All galaxy processes completed."
