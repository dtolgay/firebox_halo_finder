import numpy as np 
import pandas as pd # type: ignore

def filter_particles(particles_df, x_halo, y_halo, z_halo, rvir_halo):
    """
    Filters particles within a specified radius from the center of a halo.
    This function computes the distance of each particle from the center of a halo
    and filters out particles that are outside the specified virial radius.
    Parameters:
    particles_df (pd.DataFrame): DataFrame containing particle positions with columns 'x', 'y', and 'z'.
    x_halo (float): X-coordinate of the halo center.
    y_halo (float): Y-coordinate of the halo center.
    z_halo (float): Z-coordinate of the halo center.
    rvir_halo (float): Virial radius of the halo.
    Returns:
    pd.DataFrame: DataFrame containing only the particles within the virial radius of the halo.
    """
    

    # Compute distances from center
    distance_to_center = np.sqrt(
        (particles_df['x'] - x_halo)**2 +
        (particles_df['y'] - y_halo)**2 +
        (particles_df['z'] - z_halo)**2
    )

    # Apply filtering condition
    filtered_particles_df = particles_df[distance_to_center < rvir_halo].copy() 
    
    return filtered_particles_df

