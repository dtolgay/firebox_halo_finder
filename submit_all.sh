#!/bin/bash

# Set base values
z=0.0
start=0
end=1000
step=100

# Loop over the desired range
for ((i=start; i<end; i+=step)); do
    next=$((i + step))
    job_file="job_z${z}_${i}_${next}.sh"

    # Generate the individual job script
    cp run_python.sh "$job_file"
    sed -i "s/__Z__/$z/" "$job_file"
    sed -i "s/__START__/$i/" "$job_file"
    sed -i "s/__END__/$next/" "$job_file"
    sed -i "s/__NAME__/z2_firebox_${i}_${next}/" "$job_file"

    # Submit the job
    sbatch "$job_file"
done
