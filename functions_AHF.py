import sys
sys.path.append('/scratch/m/murray/dtolgay')
from tools import functions # type: ignore

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


####### Filter and rotate the galaxy #######


def finding_the_average_velocity_vector(vx:np.ndarray, 
                                        vy:np.ndarray,
                                        vz:np.ndarray)->list:
    
    '''This function finds the average velocity vector for given particles.
    
    Arguments:
    ----------
    vx : np.ndarray
    Velocities of the particles along the X direction

    vy : np.ndarray
    Velocities of the particles along the Y direction
    
    vz : np.ndarray
    Velocities of the particles along the Z direction    

    Returns:
    ----------
    [vx_average, vy_average, vz_average]: list
    Average velocity vectors for X, Y, and Z direction. The units are the same
    with the units of the inputted velocity arrays.
    
    '''
    vx_average = np.sum(vx)/len(vx)
    vy_average = np.sum(vy)/len(vy)
    vz_average = np.sum(vz)/len(vz)
    
    return [vx_average, vy_average, vz_average]

def shifting_velocities(vx:np.ndarray,
                        vy:np.ndarray,
                        vz:np.ndarray,
                        v_average:list)->(np.ndarray,
                                          np.ndarray,
                                          np.ndarray):
    
    '''This function shifts the velocities according to the average velocities list found above.
    
    Arguments:
    ----------
    vx : np.ndarray
    Velocities of the particles along the X direction

    vy : np.ndarray
    Velocities of the particles along the Y direction
    
    vz : np.ndarray
    Velocities of the particles along the Z direction    

    v_average : list
    Average velocities of the particles in X, Y and Z direction    

    Returns:
    ----------
    vx : np.ndarray
    Velocities of the particles along the X direction after shifting.
    Units are the same with the inputted velocity units.

    vy : np.ndarray
    Velocities of the particles along the Y direction after shifting.
    Units are the same with the inputted velocity units.
    
    vz : np.ndarray
    Velocities of the particles along the Z direction after shifting.
    Units are the same with the inputted velocity units.
    
    '''
    
    vx -= v_average[0] 
    vy -= v_average[1] 
    vz -= v_average[2] 
    
    return (vx, vy, vz)

def rotate_galaxy(gas_particles_df, star_particles_df):
    
    L_gas = functions.net_angular_momentum(
        mass=gas_particles_df["mass"], 
        rx=gas_particles_df["x"],
        ry=gas_particles_df["y"], 
        rz=gas_particles_df["z"], 
        vx=gas_particles_df["vx"], 
        vy=gas_particles_df["vy"], 
        vz=gas_particles_df["vz"]
    ) 
    # lx, ly, and lz are the indices 0, 1, and 2 respectively the net angular momentum of gas particles
    # L unit: [1e10 M☉ kpc km / sec]


        # Finding the angles between coordinate axis and net angular momentum
    theta, phi = functions.finding_the_angles_between_current_coordinate_system_and_net_angular_momentum(L=L_gas)
    # theta [radian]
    # phi   [radian]


        # Rotating the coordinate system such that z axis of the net angular momentum coincides with the positive z axis of the new coordinate system
    x_star, y_star, z_star = functions.rotating_coordinate_system_along_net_angular_momentum(
        theta=theta, 
        phi=phi, 
        vectorx=star_particles_df["x"], 
        vectory=star_particles_df["y"], 
        vectorz=star_particles_df["z"]
    )


    vx_star, vy_star, vz_star = functions.rotating_coordinate_system_along_net_angular_momentum(
        theta=theta, 
        phi=phi, 
        vectorx=star_particles_df["vx"], 
        vectory=star_particles_df["vy"], 
        vectorz=star_particles_df["vz"]
    )


    x_gas, y_gas, z_gas = functions.rotating_coordinate_system_along_net_angular_momentum(
        theta=theta, 
        phi=phi,
        vectorx=gas_particles_df["x"],
        vectory=gas_particles_df["y"],
        vectorz=gas_particles_df["z"]
    )

    vx_gas, vy_gas, vz_gas = functions.rotating_coordinate_system_along_net_angular_momentum(
        theta=theta, 
        phi=phi,
        vectorx=gas_particles_df["vx"],
        vectory=gas_particles_df["vy"],
        vectorz=gas_particles_df["vz"]
    )


    # Editing the dataframes 
    gas_particles_df["x"] = x_gas
    gas_particles_df["y"] = y_gas
    gas_particles_df["z"] = z_gas
    gas_particles_df["vx"] = vx_gas
    gas_particles_df["vy"] = vy_gas
    gas_particles_df["vz"] = vz_gas

    star_particles_df["x"] = x_star
    star_particles_df["y"] = y_star
    star_particles_df["z"] = z_star
    star_particles_df["vx"] = vx_star
    star_particles_df["vy"] = vy_star
    star_particles_df["vz"] = vz_star    
    
    return gas_particles_df, star_particles_df


def process_and_rotate_galaxy(gas_particles_df, star_particles_df):

    v_average = finding_the_average_velocity_vector(
        vx=gas_particles_df["vx"],
        vy=gas_particles_df["vy"],
        vz=gas_particles_df["vz"]
    )


    # The average velocity is subtracted from the velocities of the particles individually
    vx_gas, vy_gas, vz_gas = shifting_velocities(
        vx=gas_particles_df["vx"], 
        vy=gas_particles_df["vy"], 
        vz=gas_particles_df["vz"], 
        v_average=v_average
    )    
    
    # The average velocity is subtracted from the velocities of the particles individually
    vx_star, vy_star, vz_star = shifting_velocities(
        vx=star_particles_df["vx"], 
        vy=star_particles_df["vy"], 
        vz=star_particles_df["vz"], 
        v_average=v_average
    )

    gas_particles_df["vx"] = vx_gas
    gas_particles_df["vy"] = vy_gas
    gas_particles_df["vz"] = vz_gas  
    
    star_particles_df["vx"] = vx_star
    star_particles_df["vy"] = vy_star
    star_particles_df["vz"] = vz_star                
    
    #############################################################################################################################################
    '''
    Find the indices of the star and gas particles that are within the 20 kpc from the center of the halo 
    '''
    # Filtering the code to increase its speed
        # To increase the speed of the code I will only consider 20 kpc radius
    R_max = 20.0

    print(f"Considering only {R_max} kpc from the center of the MMH.")
    print(f"Before: len(gas_particles_df): {len(gas_particles_df)} --- len(star_particles_df): {len(star_particles_df)}")

        # Finding the distance of gas particles from the center of the MMH
    R_gas   = np.sqrt(np.power(gas_particles_df["x"],2) + np.power(gas_particles_df["y"],2) + np.power(gas_particles_df["z"],2))
    R_star  = np.sqrt(np.power(star_particles_df["x"],2) + np.power(star_particles_df["y"],2) + np.power(star_particles_df["z"],2))

        # Determining the indices that satisfy R_gas < R_max
    R_gas_smaller_than_Rmax_indices = np.where(R_gas < R_max)[0]

        # Determining the indices that satisfy R_star < R_max
    R_star_smaller_than_Rmax_indices = np.where(R_star < R_max)[0]

    # """
    # Filter the gas and star particles such that at the end you will end up with only the particles within the 20 kpc from the center
    # of the halo
    # """

    gas_particles_df = gas_particles_df.iloc[R_gas_smaller_than_Rmax_indices].reset_index(drop=True)
    star_particles_df = star_particles_df.iloc[R_star_smaller_than_Rmax_indices].reset_index(drop=True)

    print(f"After: len(gas_particles_df): {len(gas_particles_df)} --- len(star_particles_df): {len(star_particles_df)}")

    # Rotate the galaxy
    print(f"Rotating galaxy")
    gas_particles_df, star_particles_df = rotate_galaxy(
        gas_particles_df.copy(), 
        star_particles_df.copy()
    )

    return  gas_particles_df, star_particles_df