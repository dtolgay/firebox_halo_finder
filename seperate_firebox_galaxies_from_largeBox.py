import sys 
sys.path.append("/scratch/m/murray/dtolgay/")
from tools import readsnap_fireboxSeperatingGalaxies # type: ignore
import numpy as np # type: ignore 
import pandas as pd  # type: ignore
import h5py
import functions_AHF

def main(redshift):

    if redshift == "0.0":
        snapshot_number = 1200     # z = 0.0
        ahf_file_name = "FB15N1024.z0.000.AHF_halos"        
    elif redshift == "0.5":
        print(f"Exiting... Currently there are no z=0.5 galaxies... {redshift}")
        sys.exit(2)                
    elif redshift == "1.0":
        snapshot_number = 554     # z = 1.0
        ahf_file_name = "FB15N1024.z1.000.AHF_halos"        
    elif redshift == "2.0":
        snapshot_number = 344     # z = 2.0
        ahf_file_name = "FB15N1024.z2.000.AHF_halos"        
    elif redshift == "3.0":
        snapshot_number = 240     # z = 3.0
        ahf_file_name = "FB15N1024.z3.000.AHF_halos"        
    else:
        print(f"Exiting... Redshift is wrong. The given galaxy type is {redshift}")
        sys.exit(2)         

    # Define the path to the snapshot
    snap_dir_file_path = "/scratch/m/murray/dtolgay/firebox/FB15N1024"
    snapshot_number = '%03d' %snapshot_number

    # Read gas and star particles 
    gas_particles  = readsnap_fireboxSeperatingGalaxies.readsnap(snap_dir_file_path, snapshot_number, 0, cosmological=1)
    star_particles = readsnap_fireboxSeperatingGalaxies.readsnap(snap_dir_file_path, snapshot_number, 4, cosmological=1)

    # TODO: Delete below 
    print("gas and star particles from FIRE read.")

    # Create dataframe for gas particles
    gas_particles_df = pd.DataFrame({
        'x': gas_particles['p'][:,0],  # kpc 
        'y': gas_particles['p'][:,1],  # kpc 
        'z': gas_particles['p'][:,2],  # kpc 
        'vx': gas_particles['v'][:,0],
        'vy': gas_particles['v'][:,1],
        'vz': gas_particles['v'][:,2],
        'mass': gas_particles['m'],
        'density': gas_particles['rho'],
        'smoothing_length': gas_particles['h'],
        'star_formation_rate': gas_particles['sfr'],
        'internal_energy': gas_particles['u'] * 1e6,  # Converted to [m^2 s^-2]
        'neutral_hydrogen_fraction': gas_particles['nh'],
        'electron_abundance': gas_particles['ne'],
        # 'metallicity': gas_particles['z'][:,0],
        'metallicity': gas_particles['z'],
        # 'He_mass_fraction': gas_particles['z'][:,1]
        # You can add other fractions as needed
    })


    # Create dataframe for star particles
    star_particles_df = pd.DataFrame({
        'x': star_particles['p'][:,0],
        'y': star_particles['p'][:,1],
        'z': star_particles['p'][:,2],
        'vx': star_particles['v'][:,0],
        'vy': star_particles['v'][:,1],
        'vz': star_particles['v'][:,2],
        # 'metallicity': star_particles['z'][:,0],
        'metallicity': star_particles['z'],
        'scale_factor': star_particles['age'],
        'mass': star_particles['m']
    })

    # TODO: Delete below 
    print("Dataframes for gas and star particles are created.")

    # Read header info
    header_info = readsnap_fireboxSeperatingGalaxies.readsnap(snap_dir_file_path, snapshot_number, 0, header_only=1)
    hubble      = float(header_info['hubble'])
    redshift_from_header    = float(header_info['redshift'])   
    time        = float(header_info['time'])   

    # Read AHF files 
    fdir = "/scratch/m/murray/dtolgay/firebox/FB15N1024/analysis/AHF/halo_new/1200"
    

    column_names = [
        "ID", "hostHalo", "numSubStruct", "Mvir", "npart", "Xc", "Yc", "Zc", "VXc", "VYc", "VZc",
        "Rvir", "Rmax", "r2", "mbp_offset", "com_offset", "Vmax", "v_esc", "sigV", "lambda", "lambdaE",
        "Lx", "Ly", "Lz", "b", "c", "Eax", "Eay", "Eaz", "Ebx", "Eby", "Ebz", "Ecx", "Ecy", "Ecz",
        "ovdens", "nbins", "fMhires", "Ekin", "Epot", "SurfP", "Phi0", "cNFW", "n_gas", "M_gas",
        "lambda_gas", "lambdaE_gas", "Lx_gas", "Ly_gas", "Lz_gas", "b_gas", "c_gas", "Eax_gas",
        "Eay_gas", "Eaz_gas", "Ebx_gas", "Eby_gas", "Ebz_gas", "Ecx_gas", "Ecy_gas", "Ecz_gas",
        "Ekin_gas", "Epot_gas", "n_star", "M_star", "lambda_star", "lambdaE_star", "Lx_star",
        "Ly_star", "Lz_star", "b_star", "c_star", "Eax_star", "Eay_star", "Eaz_star", "Ebx_star",
        "Eby_star", "Ebz_star", "Ecx_star", "Ecy_star", "Ecz_star", "Ekin_star", "Epot_star"
    ]

    halos = pd.DataFrame(
        np.loadtxt(f"{fdir}/{ahf_file_name}"),
        columns=column_names
    )

    # Reorder the dataframe according to the Mvir 
    halos = halos.sort_values(by="Mvir", ascending=False)    


    ############### Create hdf5 files  ###############
    # Initiate an empty array to store the halos
    halos_of_interest = []

    # Find the first 100 halos 
    for row, halo in halos.iloc[0:2].iterrows():

        x_center = halo['Xc'] # kpc/h
        y_center = halo['Yc'] # kpc/h
        z_center = halo['Zc'] # kpc/h
        rvir = halo['Rvir'] # kpc/h
        mvir = halo["Mvir"]
        
        # The positions and mass of the Halo is in comoving units. To make everything on physical units the below code is written
        x_center = x_center * time/hubble     # [kpc]
        y_center = y_center * time/hubble     # [kpc]
        z_center = z_center * time/hubble     # [kpc]  
        rvir = rvir * time/hubble     # [kpc]      
        mvir = mvir * 1/hubble # [Msolar]
        
        halos_of_interest.append(
            {
                "x": x_center,
                "y": y_center,
                "z": z_center,
                "rvir": rvir,
                "mvir": mvir,
            }
        )
        
        # TODO: Change it to rvir 
        rmax = 50 # kpc 
        filtered_gas_particles = functions_AHF.filter_particles(particles_df = gas_particles_df, x_halo = x_center, y_halo = y_center, z_halo = z_center, rmax = rmax)
        filtered_star_particles = functions_AHF.filter_particles(particles_df = star_particles_df, x_halo = x_center, y_halo = y_center, z_halo = z_center, rmax = rmax)

        # Change the origin of the particles 
        filtered_gas_particles = functions_AHF.change_origin(
            particles = filtered_gas_particles.copy(),
            x_halo_center = x_center, 
            y_halo_center = y_center, 
            z_halo_center = z_center,
        )
        
        filtered_star_particles = functions_AHF.change_origin(
            particles = filtered_star_particles.copy(), 
            x_halo_center = x_center, 
            y_halo_center = y_center,
            z_halo_center = z_center, 
        )
        
        # Reconstruct the dictionary to write the file 
        reconstructed_gas_particles = {
            'p': filtered_gas_particles[['x', 'y', 'z']].to_numpy(),  # shape (N, 3)
            'v': filtered_gas_particles[['vx', 'vy', 'vz']].to_numpy(),  # shape (N, 3)
            'm': filtered_gas_particles['mass'].to_numpy(),  # shape (N,)
            'rho': filtered_gas_particles['density'].to_numpy(),
            'h': filtered_gas_particles['smoothing_length'].to_numpy(),
            'sfr': filtered_gas_particles['star_formation_rate'].to_numpy(),
            'u': (filtered_gas_particles['internal_energy'] / 1e6).to_numpy(),  # reverse conversion
            'nh': filtered_gas_particles['neutral_hydrogen_fraction'].to_numpy(),
            'ne': filtered_gas_particles['electron_abundance'].to_numpy(),
            # 'z': filtered_gas_particles[['metallicity', 'He_mass_fraction']].to_numpy()  # shape (N, 2)
            'z': filtered_gas_particles['metallicity'].to_numpy()  
        }
        
        reconstructed_star_particles = {
            'p': filtered_star_particles[['x', 'y', 'z']].to_numpy(),         # shape (N, 3)
            'v': filtered_star_particles[['vx', 'vy', 'vz']].to_numpy(),      # shape (N, 3)
            'z': filtered_star_particles[['metallicity']].to_numpy(),         # shape (N, 1) — or (N,) depending on original
            'age': filtered_star_particles['scale_factor'].to_numpy(),        # shape (N,)
            'm': filtered_star_particles['mass'].to_numpy()                   # shape (N,)
        }    
        
        # Save galaxies into a hdf5 file
        file_name = f"gal_{row}.hdf5"
        with h5py.File(f'z{redshift}/{file_name}', 'w') as f: 
            gas_group = f.create_group('gas')
            for key, value in reconstructed_gas_particles.items():
                gas_group.create_dataset(key, data=value)

            stars_group = f.create_group('star')
            for key, value in reconstructed_star_particles.items():
                stars_group.create_dataset(key, data=value)
        
        print(f"{file_name} written!")
                

    # TODO: I am doing nothing for now with the halos_of_interest.
    halos_of_interest = pd.DataFrame(halos_of_interest)
        
    return None 


if __name__ == "__main__":
    redshift = sys.argv[1]
    main(redshift=redshift)