import numpy as np 
import pandas as pd # type: ignore

def filter_particles(particles_df, x_halo, y_halo, z_halo, rmax):
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
    filtered_particles_df = particles_df[distance_to_center < rmax].copy() 
    
    return filtered_particles_df

def change_origin(particles, x_halo_center, y_halo_center, z_halo_center):

    particles_new = particles.copy()

    particles_new['x'] = particles_new['x'] - x_halo_center 
    particles_new['y'] = particles_new['y'] - y_halo_center 
    particles_new['z'] = particles_new['z'] - z_halo_center 

    return particles_new

def create_df_from_read_hdf5(particles, particle_type):
    """
    Create a pandas DataFrame from HDF5 particle data based on the particle type.

    Parameters:
    -----------
    particles : dict
        A dictionary containing particle data arrays. The keys should correspond to the 
        particle properties such as position, velocity, mass, etc.
    particle_type : str
        The type of particles to process. Must be either 'gas' or 'star'.

    Returns:
    --------
    pandas.DataFrame
        A DataFrame containing the particle data with appropriate columns based on the 
        particle type.

    Raises:
    -------
    ValueError
        If an invalid particle type is provided.
    """

    if particle_type == "gas":
        particles_df = pd.DataFrame({
            'x': particles['p'][:,0],                   # kpc 
            'y': particles['p'][:,1],                   # kpc 
            'z': particles['p'][:,2],                   # kpc 
            'vx': particles['v'][:,0],                  # km/s
            'vy': particles['v'][:,1],                  # km/s
            'vz': particles['v'][:,2],                  # km/s
            'mass': particles['m'],                     # 1e10 Msun
            'density': particles['rho'],                # 1e10 Msun / kpc^3
            'smoothing_length': particles['h'],         # kpc
            'star_formation_rate': particles['sfr'],    # Msun / yr
            'internal_energy': particles['u'] * 1e6,  
            'neutral_hydrogen_fraction': particles['nh'],  
            'electron_abundance': particles['ne'],
            'metallicity': particles['z'][:,0],         # mass fraction
            'He_mass_fraction': particles['z'][:,1]
        })        
    elif particle_type ==  "star":
        particles_df = pd.DataFrame({
            'x': particles['p'][:,0],
            'y': particles['p'][:,1],
            'z': particles['p'][:,2],
            'vx': particles['v'][:,0],
            'vy': particles['v'][:,1],
            'vz': particles['v'][:,2],
            'metallicity': particles['z'][:,0],
            'scale_factor': particles['age'],
            'mass': particles['m']
        })
    else:
        raise ValueError("Invalid particle type. Choose between 'gas' and 'star'.")

    return particles_df 

