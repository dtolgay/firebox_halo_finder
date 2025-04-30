import sys 
sys.path.append("/scratch/m/murray/dtolgay/")
from tools import readsnap_fireboxSeperatingGalaxies # type: ignore
from tools.readsnap import check_if_filename_exists, load_gadget_format_binary_header, load_gadget_format_binary_particledat # type: ignore

import numpy as np # type: ignore 
import pandas as pd  # type: ignore
import h5py # type: ignore
import functions_AHF
import traceback

def main(redshift, start_halo, stop_halo):

    if redshift == "0.0":
        snapshot_number = 1200     # z = 0.0
        ahf_file_name = "1200/FB15N1024.z0.000.AHF_halos"        
    elif redshift == "0.5":
        print(f"Exiting... Currently there are no z=0.5 galaxies... {redshift}")
        sys.exit(2)                
    elif redshift == "1.0":
        snapshot_number = 554     # z = 1.0
        ahf_file_name = "554/FB15N1024.z1.000.AHF_halos"        
    elif redshift == "2.0":
        snapshot_number = 344     # z = 2.0
        ahf_file_name = "344/FB15N1024.z2.000.AHF_halos"        
    elif redshift == "3.0":
        snapshot_number = 240     # z = 3.0
        ahf_file_name = "240/FB15N1024.z3.000.AHF_halos"        
    else:
        print(f"Exiting... Redshift is wrong. The given galaxy type is {redshift}")
        sys.exit(2)       

    #### Read header ####
    # Define the path to the snapshot
    snap_dir_file_path = "/scratch/m/murray/dtolgay/firebox/FB15N1024"    
    header_info = readsnap_only_one_file_at_a_time(sdir=snap_dir_file_path, snum=snapshot_number, ptype=0, which_snapshot_to_read=0, header_only=1)
    hubble      = float(header_info['hubble'])
    redshift_from_header    = float(header_info['redshift'])   
    time        = float(header_info['time'])       
    NumFilesPerSnapshot = header_info['NumFilesPerSnapshot']


    ##### Read halos #####
    # Read AHF files 
    fdir = "/scratch/m/murray/dtolgay/firebox/FB15N1024/analysis/AHF/halo_new"
    
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

    halo_fdir = f"{fdir}/{ahf_file_name}"
    halos = pd.DataFrame(
        np.loadtxt(halo_fdir),
        columns=column_names
    )
    print(f"{halo_fdir} is read!")

    # Sort the halos according to the Mvir in descending order    
    halos = halos.sort_values(by="Mvir", ascending=False).reset_index(drop=True) 
    halos['new_id_after_sorting'] = halos.index
    
    # Write all halos into a file 
    halos.to_csv(f"z{redshift}/halos_used.csv", index=False)
    print(f"All halos are written to z{redshift}/halos_used.csv")

    ### For loop to seperete galaxies from the large box ###
    total_number_of_gas_particles = 0
    halos_of_interest = []
    for row, halo in halos.iloc[start_halo:stop_halo].iterrows(): # TODO:
        print("\n\n\n")
        print(f" ----------- Halo: {row} ----------- ")

        # if hostHalo is zero, that means it is the top level halo
        if (int(halo['hostHalo']) == -1): 
            
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
            
            halos_of_interest.append(halo) # Append halo for later use 

            ### Read the snapshot and find the gas particles inside the halos ###
            try:
                # Reading the particles in each snapshot seperately. 
                # Initiate the filtered gas and star particles dataframes
                filtered_gas_particles = pd.DataFrame()
                filtered_star_particles = pd.DataFrame()
                for i in range(NumFilesPerSnapshot):  
                    print(f" ----------- which snapshot: {i} ----------- ")
                    # Read gas and star particles
                    gas_particles  = readsnap_only_one_file_at_a_time(sdir=snap_dir_file_path, snum=snapshot_number, ptype=0, which_snapshot_to_read=i, cosmological=1)
                    star_particles = readsnap_only_one_file_at_a_time(sdir=snap_dir_file_path, snum=snapshot_number, ptype=4, which_snapshot_to_read=i, cosmological=1)
                    total_number_of_gas_particles += len(gas_particles['p'][:,0])

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
                        'metallicity': gas_particles['z'][:,0],
                        'He_mass_fraction': gas_particles['z'][:,1],
                        'C_mass_fraction': gas_particles['z'][:,2],
                        'N_mass_fraction': gas_particles['z'][:,3],
                        'O_mass_fraction': gas_particles['z'][:,4],
                        'Ne_mass_fraction': gas_particles['z'][:,5],
                        'Mg_mass_fraction': gas_particles['z'][:,6],
                        'Si_mass_fraction': gas_particles['z'][:,7],
                        'S_mass_fraction': gas_particles['z'][:,8],
                        'Ca_mass_fraction': gas_particles['z'][:,9],
                        'Fe_mass_fraction': gas_particles['z'][:,10],
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
                        'metallicity': star_particles['z'][:,0],
                        'He_mass_fraction': star_particles['z'][:,1],
                        'C_mass_fraction': star_particles['z'][:,2],
                        'N_mass_fraction': star_particles['z'][:,3],
                        'O_mass_fraction': star_particles['z'][:,4],
                        'Ne_mass_fraction': star_particles['z'][:,5],
                        'Mg_mass_fraction': star_particles['z'][:,6],
                        'Si_mass_fraction': star_particles['z'][:,7],
                        'S_mass_fraction': star_particles['z'][:,8],
                        'Ca_mass_fraction': star_particles['z'][:,9],
                        'Fe_mass_fraction': star_particles['z'][:,10],                        
                        'scale_factor': star_particles['age'],
                        'mass': star_particles['m']
                    })                


                    # Find the gas particles inside the halos 
                    rmax = 50 # kpc 
                    filtered_gas_particles_for_this_snapshot = functions_AHF.filter_particles(particles_df = gas_particles_df, x_halo = x_center, y_halo = y_center, z_halo = z_center, rmax = rmax)
                    filtered_star_particles_for_this_snapshot = functions_AHF.filter_particles(particles_df = star_particles_df, x_halo = x_center, y_halo = y_center, z_halo = z_center, rmax = rmax)                


                    # Change the origin of the particles 
                    filtered_gas_particles_for_this_snapshot = functions_AHF.change_origin(
                        particles = filtered_gas_particles_for_this_snapshot.copy(),
                        x_halo_center = x_center, 
                        y_halo_center = y_center, 
                        z_halo_center = z_center,
                    )

                    filtered_star_particles_for_this_snapshot = functions_AHF.change_origin(
                        particles = filtered_star_particles_for_this_snapshot.copy(), 
                        x_halo_center = x_center, 
                        y_halo_center = y_center,
                        z_halo_center = z_center, 
                    )


                    # Append filteres gas and star particles to their respective dataframes
                    filtered_gas_particles = pd.concat([filtered_gas_particles, filtered_gas_particles_for_this_snapshot])
                    filtered_star_particles = pd.concat([filtered_star_particles, filtered_star_particles_for_this_snapshot])

                # Recreate the dictionary datas tructure to write gas and star particles to hdf5 file 
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
                    'z': filtered_gas_particles[
                        [
                            'metallicity', 
                            'He_mass_fraction', 
                            'C_mass_fraction',
                            'N_mass_fraction',
                            'O_mass_fraction',
                            'Ne_mass_fraction',
                            'Mg_mass_fraction',
                            'Si_mass_fraction',
                            'S_mass_fraction',
                            'Ca_mass_fraction',
                            'Fe_mass_fraction' 
                        ]].to_numpy()                                                 # shape (N, 11)
                }

                reconstructed_star_particles = {
                    'p': filtered_star_particles[['x', 'y', 'z']].to_numpy(),         # shape (N, 3)
                    'v': filtered_star_particles[['vx', 'vy', 'vz']].to_numpy(),      # shape (N, 3)
                    'z': filtered_star_particles[
                        [
                            'metallicity', 
                            'He_mass_fraction', 
                            'C_mass_fraction',
                            'N_mass_fraction',
                            'O_mass_fraction',
                            'Ne_mass_fraction',
                            'Mg_mass_fraction',
                            'Si_mass_fraction',
                            'S_mass_fraction',
                            'Ca_mass_fraction',
                            'Fe_mass_fraction' 
                        ]].to_numpy(),                                                # shape (N, 11)                    
                    'age': filtered_star_particles['scale_factor'].to_numpy(),        # shape (N,)
                    'm': filtered_star_particles['mass'].to_numpy()                   # shape (N,)
                }

                # Save galaxies into a hdf5 file
                file_name = f"gal_{row}.hdf5"
                with h5py.File(f'z{redshift}/{file_name}', 'w') as f:
                    # Gas 
                    gas_group = f.create_group('gas')
                    for key, value in reconstructed_gas_particles.items():
                        gas_group.create_dataset(key, data=value)

                    # Stars
                    stars_group = f.create_group('star')
                    for key, value in reconstructed_star_particles.items():
                        stars_group.create_dataset(key, data=value)

                    # header 
                    header_group = f.create_group('header')
                    for key, value in header_info.items():
                        header_group.create_dataset(key, data=value)

                print(f"{file_name} written!")                            

            except Exception as e:
                print("Exception occurred:")
                print(f"Type: {type(e)}")
                print(f"Message: {e}")
                traceback.print_exc()


            # Write halos to a file 
            halos_of_interest_df = pd.DataFrame(halos_of_interest)
            halos_of_interest_df.to_csv(f"z{redshift}/halos_used_{start_halo}_{stop_halo}.csv", index=False)
            

        else: 
            print("This is not the top level halo.")
            print(halo)
            continue


        print(f"Total number of gas particles: {total_number_of_gas_particles}")


    return None 


# Functions defined below 


def readsnap_only_one_file_at_a_time(sdir,snum,ptype,
    which_snapshot_to_read,
    snapshot_name='snapshot',
    extension='.hdf5',
    h0=0,cosmological=0,skip_bh=0,four_char=0,
    header_only=0,loud=0):
    '''
    This is a sub-routine designed to copy a GIZMO snapshot portion - specifically
    all the data corresponding to particles of a given type - into active memory in 
    a parent structure for use in python. The routine automatically handles multi-part 
    snapshot files for you (concatenating), and works with python2.x and python3.x, 
    and both GIZMO hdf5 and un-formatted binary outputs. -- dtolgay: This is not longer true. It only reads one snapshot file at a time.

    Syntax:
      P = readsnap(sdir,snum,ptype,....)
      
      Here "P" is a structure which contains the data. The snapshot file[s] are opened,
      the data fully copied out, and the file closed. This attempts to copy 
      all the common data types, all together, into "P". Three things to note: 
      (1) the fields in P (visible by typing P.keys()) are not given the same names 
      as those in the raw snapshot, but 'shorthand' names, for convenience. you should
      look at the keys and be sure you know which associate with which files.
      (2) because of the full-copy approach into new-named fields, this will not 
      handle arbitrary new data types (it is impossible to handle these in full 
      generality with un-formatted binary, you need to code the file-order and byte 
      numbers for each new data structure). if you add new fields (with e.g. 
      additional physics modules) beyond what this code looks for, you need to 
      add code here, or use the more general 'load_from_snapshot.py' routine.
      (3) also because of the full copy strategy, this routine is much more expensive 
      in time and memory compared to 'load_from_snapshot.py'. use that if you want 
      a light-weight, more flexible reading option, and have HDF5 outputs.

      For example, after calling, the 'Coordinates' field from the snapshot is 
      accessible from the new structure P by calling P['p']. 'Velocities' as P['v'],
      'Masses' as P['m']. Fields specific to gas include 'Density' as P['rho'], 
      'InternalEnergy' as P['u'], and more.

      More details and examples are given in the GIZMO user guide.

    Arguments:               
      sdir: parent directory (string) of the snapshot file or immediate snapshot sub-directory 
            if it is a multi-part file.
            
      snum: number (int) of the snapshot. e.g. snapshot_001.hdf5 is '1'
            Note for multi-part files, this is just the number of the 'set', i.e. 
            if you have snapshot_001.N.hdf5, set this to '1', not 'N' or '1.N'

      ptype: element type (int) = 0[gas],1,2,3,4,5[meaning depends on simulation, see
             user guide for details]. if your chosen 'value' is in the file header, 
             this will be ignored
      
    Optional:
      header_only: the structure "P" will return the file header, instead of the 
        particle data. you can see the data in the header then by simply typing 
        P.keys() -- this contains data like the time of the snapshot, as 
        P['Time']. With this specific routine the header information is only 
        saved if you choose this option. Default 0/False, turn on by setting to 
        1 or True.
      
      cosmological: default 0/False: turn on (set to 1/True) to convert cosmological 
        co-moving units to physical units. will specifically convert Coordinates, 
        Masses, Velocities, Densities, Smoothing Lengths, and Times/Ages. If this 
        is on, you do not need to set 'h0' (this will force it to be set -also-), 
        but it does no harm to set it as well.

      h0: default 0/False: turn on (set to 1/True) for the code to use the value 
        of the hubble constant (h = H0/100 km/s/Mpc) saved in the snapshot to convert 
        from code units. Recall, units of time, length, and mass in the code are in 
        h^-1. So this on means your units are physical, with no "h" in them. 
        Will specifically convert Coordinates, Masses, Densities, Smoothing Lengths, 
        and Times/Ages. 

      skip_bh: default 0/False: turn on (set to 1/True) to skip black hole-specific 
        fields for particles of type 5 (use if your snapshot contains elements of type 5, 
        but the black hole physics modules were not actually used; otherwise you will 
        get an error).
      
      four_char: default numbering is that snapshots with numbers below 1000 have 
        three-digit numbers. if they were numbered with four digits (e.g. snapshot_0001), 
        set this to 1 or True (default 0/False)
        
      snapshot_name: default 'snapshot': the code will automatically try a number of 
        common snapshot and snapshot-directory prefixes. but it can't guess all of them, 
        especially if you use an unusual naming convention, e.g. naming your snapshots 
        'xyzBearsBeetsBattleStarGalactica_001.hdf5'. In that case set this to the 
        snapshot name prefix (e.g. 'xyzBearsBeetsBattleStarGalactica')
        
      extension: default 'hdf5': again like 'snapshot' set if you use a non-standard 
        extension (it checks multiply options like 'h5' and 'hdf5' and 'bin'). but 
        remember the file must actually be hdf5 format!

      loud: print additional checks as it reads, useful for debugging, 
        set to 1 or True if desired (default 0/False)
    


    '''

    if (ptype<0): return {'k':-1};
    if (ptype>5): return {'k':-1};

    fname_for_snaphot_zero,fname_base,fname_ext = check_if_filename_exists(sdir,snum,\
        snapshot_name=snapshot_name,extension=extension,four_char=four_char)
    
    fname = f"{fname_base}.{which_snapshot_to_read}{fname_ext}"
    
    if(fname=='NULL'): return {'k':-1}
    if(loud==1): print('loading file : '+fname)

    ## open file and parse its header information
    nL = 0 # initial particle point to start at 
    if(fname_ext=='.hdf5'):
        file = h5py.File(fname,'r') # Open hdf5 snapshot file
        header_topdict = file["Header"] # Load header dictionary (to parse below)
        header_toparse = header_topdict.attrs
    else:
        file = open(fname) # Open binary snapshot file
        header_toparse = load_gadget_format_binary_header(file)

    npart = header_toparse["NumPart_ThisFile"]
    massarr = header_toparse["MassTable"]
    time = header_toparse["Time"]
    redshift = header_toparse["Redshift"]
    flag_sfr = header_toparse["Flag_Sfr"]
    flag_feedbacktp = header_toparse["Flag_Feedback"]
    npartTotal = header_toparse["NumPart_Total"]
    flag_cooling = header_toparse["Flag_Cooling"]
    numfiles = header_toparse["NumFilesPerSnapshot"]
    boxsize = header_toparse["BoxSize"]
    hubble = header_toparse["HubbleParam"]
    flag_stellarage = header_toparse["Flag_StellarAge"]
    flag_metals = header_toparse["Flag_Metals"]
    print("npart_file: ",npart)
    print("npart_total:",npartTotal)

    hinv=1.
    if (h0==1):
        hinv=1./hubble
    ascale=1.
    if (cosmological==1):
        ascale=time
        hinv=1./hubble
    if (cosmological==0): 
        time*=hinv
    
    
    boxsize*=hinv*ascale
    if (npartTotal[ptype]<=0): file.close(); return {'k':-1};
    if (header_only==1): file.close(); return {'k':0,'time':time,
        'boxsize':boxsize,'hubble':hubble,'npart':npart,'npartTotal':npartTotal, 'redshift':redshift, 'NumFilesPerSnapshot':numfiles};
    # Update: redshift is added by dtolgay

    # initialize variables to be read
    pos=np.zeros([npart[ptype],3],dtype=np.float64) 
    vel=np.copy(pos)
    ids=np.zeros([npart[ptype]],dtype=int)
    mass=np.zeros([npart[ptype]],dtype=np.float64)
    if (ptype==0):
        ugas=np.copy(mass)
        rho=np.copy(mass)
        hsml=np.copy(mass) 
        #if (flag_cooling>0): 
        nume=np.copy(mass)
        numh=np.copy(mass)
        #if (flag_sfr>0): 
        sfr=np.copy(mass)
        metal=np.copy(mass)
    if (ptype == 0 or ptype == 4) and (flag_metals > 0):
        metal=np.zeros([npart[ptype],flag_metals],dtype=np.float32)
    if (ptype == 4) and (flag_sfr > 0) and (flag_stellarage > 0):
        stellage=np.copy(mass)
    if (ptype == 5) and (skip_bh == 0):
        bhmass=np.copy(mass)
        bhmdot=np.copy(mass)
    
    # open the file            
    if (fname_ext=='.hdf5'):
        input_struct = file
        npart = file["Header"].attrs["NumPart_ThisFile"]
        bname = "PartType"+str(ptype)+"/"
    else:
        npart = header_toparse['NumPart_ThisFile']
        input_struct = load_gadget_format_binary_particledat(file, header_toparse, ptype, skip_bh=skip_bh)
        bname = ''
            
    
    # now do the actual reading
    if(npart[ptype]>0):
        nR=nL + npart[ptype]
        pos[nL:nR,:]=input_struct[bname+"Coordinates"]
        vel[nL:nR,:]=input_struct[bname+"Velocities"]
        ids[nL:nR]=input_struct[bname+"ParticleIDs"]
        mass[nL:nR]=massarr[ptype]
        if (massarr[ptype] <= 0.):
            mass[nL:nR]=input_struct[bname+"Masses"]
        if (ptype==0):
            ugas[nL:nR]=input_struct[bname+"InternalEnergy"]
            rho[nL:nR]=input_struct[bname+"Density"]
            hsml[nL:nR]=input_struct[bname+"SmoothingLength"]
            if (flag_cooling > 0): 
                nume[nL:nR]=input_struct[bname+"ElectronAbundance"]
                numh[nL:nR]=input_struct[bname+"NeutralHydrogenAbundance"]
            if (flag_sfr > 0):
                sfr[nL:nR]=input_struct[bname+"StarFormationRate"]
        if (ptype == 0 or ptype == 4) and (flag_metals > 0):
            metal_t=input_struct[bname+"Metallicity"]
            if (flag_metals > 1):
                if (metal_t.shape[0] != npart[ptype]): 
                    metal_t=np.transpose(metal_t)
            else:
                metal_t=np.reshape(np.array(metal_t),(np.array(metal_t).size,1))
            metal[nL:nR,:]=metal_t
        if (ptype == 4) and (flag_sfr > 0) and (flag_stellarage > 0):
            stellage[nL:nR]=input_struct[bname+"StellarFormationTime"]
        if (ptype == 5) and (skip_bh == 0):
            bhmass[nL:nR]=input_struct[bname+"BH_Mass"]
            bhmdot[nL:nR]=input_struct[bname+"BH_Mdot"]
        nL = nR # sets it for the next iteration	

	## correct to same ID as original gas particle for new stars, if bit-flip applied
    if ((np.min(ids)<0) | (np.max(ids)>1.e9)):
        bad = (ids < 0) | (ids > 1.e9)
        ids[bad] += (int(1) << 31)

    # do the cosmological conversions on final vectors as needed
    pos *= hinv*ascale # snapshot units are comoving
    mass *= hinv
    vel *= np.sqrt(ascale) # remember gizmo's (and gadget's) weird velocity units!
    if (ptype == 0):
        rho *= (hinv/((ascale*hinv)**3))
        hsml *= hinv*ascale
    if (ptype == 4) and (flag_sfr > 0) and (flag_stellarage > 0) and (cosmological == 0):
        stellage *= hinv
    if (ptype == 5) and (skip_bh == 0):
        bhmass *= hinv

    file.close();
    if (ptype == 0):
        return {'k':1,'p':pos,'v':vel,'m':mass,'id':ids,'u':ugas,'rho':rho,'h':hsml,'ne':nume,'nh':numh,'sfr':sfr,'z':metal};
    if (ptype == 4):
        return {'k':1,'p':pos,'v':vel,'m':mass,'id':ids,'z':metal,'age':stellage}
    if (ptype == 5) and (skip_bh == 0):
        return {'k':1,'p':pos,'v':vel,'m':mass,'id':ids,'mbh':bhmass,'mdot':bhmdot}
    return {'k':1,'p':pos,'v':vel,'m':mass,'id':ids}




if __name__ == "__main__":
    redshift = sys.argv[1]
    start_halo = int(sys.argv[2])
    stop_halo = int(sys.argv[3])
    main(redshift=redshift, start_halo=start_halo, stop_halo=stop_halo)