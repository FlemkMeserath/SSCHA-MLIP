import numpy as np
import sys,os
import time
import spglib as spg
import math
from ase.build import bulk

def read_dyn0(MATDYN):
    dyn0_file = open(MATDYN + str(0),"r")
    lines_to_read = dyn0_file.readlines()
    mesh = lines_to_read[0].split()
    nmats = lines_to_read[1].split()[0]
    qmesh = (int(mesh[0]),int(mesh[1]),int(mesh[2]))
    return qmesh,int(nmats)



def make_input(MATDYN):
   
    cell_alat, numbers, names, mass, structure, cell, for_cons = read_dynamical_matrix(MATDYN + "1")

    supercell = [1,1,1]
    for i in range(3):
        supercell[i] = int(round( 3.57790149085*2 / math.sqrt(cell[i][0]*cell[i][0] + cell[i][1]*cell[i][1] + cell[i][2]*cell[i][2])))

    ##########################################################
    #   PRINTING KPOINT SAMP                                 #
    ##########################################################
    testing = bulk('Au').cell
    rep_latt_lenght = [0,0,0]
    
    for i in range(0,3):
        for j in range(0,3):
            testing[i][j] = cell[i][j]*supercell[i]
    
    
    Rcell = testing.reciprocal()
    for i in range(0,3):
        SwapV = np.sqrt( Rcell[i][0]*Rcell[i][0] + Rcell[i][1]*Rcell[i][1] + Rcell[i][2]*Rcell[i][2] )
        rep_latt_lenght[i] = SwapV
    
    header = "&control\n\
        calculation        = 'scf'\n\
        restart_mode       = 'from_scratch'\n\
        prefix             = 'scf'\n\
        tstress            = .true. \n\
        tprnfor            = .true.\n\
        pseudo_dir         = './Pseudo'\n\
        outdir             = './tmp/'\n\
        verbosity          = 'high'\n\
    &end\n\
    \n\
    &system\n\
        ibrav              = 0\n\
        nat                = 216\n\
        ntyp               = 1\n\
        ecutwfc            = 70.0\n\
        ecutrho            = 700\n\
        degauss            = 0.01\n\
        occupations        = 'smearing'\n\
        smearing           = 'mp'\n\
     &end\n\
    \n\
    &electrons\n\
        mixing_beta        = 0.2\n\
        conv_thr           = 1.0d-9\n\
     &end\n\
    &ions\n\
     &end\n\
    \n\
     ATOMIC_SPECIES\n\
      Li  6.941   Li.pbe-sl-rrkjus_psl.1.0.0.UPF\n\
    \n\
    K_POINTS (automatic)\n\
      3 3 3   1  1  1\n\
    \n"

    return header

########################################################




def read_poscar():
    cell = []
    atomic_positions = []
    poscar_file = open(str("POSCAR"),"r")
    poscar_file.readline()
    poscar_file.readline()

    line = poscar_file.readline()
    cell.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
    line = poscar_file.readline()
    cell.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
    line = poscar_file.readline()
    cell.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
    atom_types = poscar_file.readline().split()
    atoms_type_numbers_degen =  poscar_file.readline().split()
    total_n_atoms = 0
    for i in range(0,len(atoms_type_numbers_degen)):
        atoms_type_numbers_degen[i] = int(atoms_type_numbers_degen[i])
        total_n_atoms = total_n_atoms + atoms_type_numbers_degen[i]
    poscar_file.readline()
    for i in range(0,total_n_atoms):
        line = poscar_file.readline()
        atomic_positions.append([ float(line.split()[0]), float(line.split()[1]), float(line.split()[2]) ])

    poscar_file.close()

    return cell, atomic_positions, atoms_type_numbers_degen, atom_types

def read_dynamical_matrix(dyn_name):

    bohr_constant = 0.5291772109
    atom_type_mass_list = []
    atom_name_list = []
    atom_coords_list = []
    atom_type_coords_list = []
    names = []
    masses = []

    for_cons = []
    dynmat = open(str(dyn_name),"r")
    file_lines = dynmat.readlines()

    primary_alat = float(file_lines[2].split()[3])
    atom_types   = int(file_lines[2].split()[0])
    number_of_atoms = int(file_lines[2].split()[1])
    basis_vec_line1 = [float(i) for i in file_lines[4].split()]
    basis_vec_line2 = [float(i) for i in file_lines[5].split()]
    basis_vec_line3 = [float(i) for i in file_lines[6].split()]

        # Convert cell parameters
    basis_vec_line1 = [i*primary_alat for i in basis_vec_line1]
    basis_vec_line2 = [i*primary_alat for i in basis_vec_line2]
    basis_vec_line3 = [i*primary_alat for i in basis_vec_line3]
    np_cell = np.array([basis_vec_line1,basis_vec_line2,basis_vec_line3])

    for i in range(7, 7 + atom_types):
        atom_type_mass_list.append((file_lines[i].replace("'",'')).split()[1:])
        atom_name_list.append((file_lines[i].replace("'",'')).split()[:2])


    for i in range(7 + atom_types, 7 + atom_types + number_of_atoms):
        atom_coords_list.append([float(n) for n in file_lines[i].split()][2:])
        atom_type_coords_list.append(file_lines[i].split()[1])
    np_atom_list = np.array(atom_coords_list)


    for i in range(0,len(atom_type_coords_list)):
        for j in range(0,len(atom_name_list)):
            if atom_type_coords_list[i] == atom_name_list[j][0]:
                names.append(atom_name_list[j][1])
                masses.append(float(atom_type_mass_list[j][1]))

    atom_type_coords_list = [int(i)-1 for i in atom_type_coords_list]
    primary_alat = primary_alat * bohr_constant
    np_atom_list = np_atom_list*primary_alat
    np_cell = np_cell*bohr_constant

    for i in range(7 + number_of_atoms + atom_types + 5, 4*number_of_atoms*number_of_atoms+7 + number_of_atoms + atom_types + 5,4):
        for j in range(1,4):
            string = file_lines[i+j].split()
            for_cons.append([ float(string[0]), float(string[1]), float(string[2]), float(string[3]), float(string[4]), float(string[5])  ])

    return primary_alat, atom_type_coords_list,names,masses,np_atom_list,np_cell,for_cons

def Generate_CFG(POPULATION,ENS_FOLDER,atom_type_table):
    Ry_to_eV = 13.6057039763
    Bohr_to_A = 0.529177210903

    cfg_file = open("MLIP_"+str(POPULATION) + ".cfg", "w")
    directory = ENS_FOLDER + str(POPULATION)
    FORCE_FILES = True
    
    energy_file = open(directory + "/energies_supercell_population" + str(POPULATION) + ".dat", "r")
    energy_lines = [l.strip() for l in energy_file.readlines()]
    energy_file.close()
    
    print(len(energy_lines),energy_lines)
    
    for j in range(1,len(energy_lines)+1):
        structure_file = open(directory + "/scf_population" + str(POPULATION) + "_" + str(j) + ".dat", "r")
        lines = [l.strip() for l in structure_file.readlines()]
        structure_file.close()
        
        try:
            force_file = open(directory + "/forces_population" + str(POPULATION) + "_" + str(j) + ".dat", "r")
            force_lines = [l.strip() for l in force_file.readlines()]
            force_file.close()
        except: FORCE_FILES = False
    
        try:
            pressure_file = open(directory + "/pressures_population" + str(POPULATION) + "_" + str(j) + ".dat", "r")
            pressure_lines = [l.strip() for l in pressure_file.readlines()]
            pressure_file.close()
        except: 
            pressure_lines = [
            "0.0000 0.0000 0.0000",
            "0.0000 0.0000 0.0000",
            "0.0000 0.0000 0.0000",
            ]
        
        SIZE = len(lines)-6
        cell = np.zeros((3, 3))
        atoms = np.zeros((SIZE,3))
        forces = np.zeros((SIZE,3))
        atm_type = [None] * SIZE
        pressures = np.zeros((3, 3))
        
        for i in range(0,3):
            cell[i, :] = [float(x) for x in lines[i+1].split()[-3:]]
            
        for i in range(0,SIZE):
            atoms[i, :] = [float(x) for x in lines[i+6].split()[-3:]]
            atm_type[i] =  lines[i+6].split()[0]
            if FORCE_FILES == True:
                forces[i,:] = [float(x) for x in force_lines[i].split()[-3:]]
                forces[i,:] = forces[i, :] * Ry_to_eV / Bohr_to_A
        
            
        for i in atm_type:
            try: atom_type_table[i]
            except: atom_type_table[i] = len(atom_type_table)    
        
        
        Volume = np.dot(cell[0],np.cross(cell[1], cell[2]))
        
        for i in range(0,3):
        
            pressures[i, :] = [float(x) for x in pressure_lines[i].split()[-3:]]    
            pressures[i, :] = pressures[i, :] * Ry_to_eV / Bohr_to_A / Bohr_to_A / Bohr_to_A* Volume
        
        cfg_file.write("BEGIN_CFG\n")
        cfg_file.write(" Size\n")
        cfg_file.write("{: >5}".format(SIZE) +"\n")
        cfg_file.write(" Supercell\n")
        for row in cell:
            cfg_file.write("    {: >13f} {: >13f} {: >13f}\n".format(*row))
        cfg_file.write(" AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n")
        for i in range(0,len(atoms)):
            if FORCE_FILES == True:
                cfg_file.write("    {: >10}".format(i+1) + "{: >5}".format(atom_type_table[atm_type[i]]) + "  {: >13f} {: >13f} {: >13f}".format(*atoms[i,:]) + "  {: >11f} {: >11f} {: >11f}\n".format(*forces[i,:]))
            else:
                cfg_file.write("    {: >10}".format(i+1) + "{: >5}".format(atom_type_table[atm_type[i]]) + "  {: >13f} {: >13f} {: >13f}".format(*atoms[i,:]) + "  {: >11f} {: >11f} {: >11f}\n".format(*[0,0,0]))
        cfg_file.write(" Energy\n")
        cfg_file.write("     {: >13f}".format(float(energy_lines[j-1])*Ry_to_eV) +"\n")
        cfg_file.write(" PlusStress:  xx          yy          zz          yz          xz          xy\n")
        cfg_file.write("    {: >12f}".format(pressures[0,0]) + "{: >12f}".format(pressures[1,1]) +  "{: >12f}".format(pressures[2,2]))
        cfg_file.write("{: >12f}".format(pressures[1,2]) + "{: >12f}".format(pressures[0,2]) +  "{: >12f}\n".format(pressures[0,1]))
        cfg_file.write(" Feature atom_type_table " + str(atom_type_table) + "\n")
        cfg_file.write(" Feature conf_number " + str(j) + "\n")  
        cfg_file.write(" Feature population " + str(POPULATION) + "\n")  
        cfg_file.write("END_CFG\n\n")
    cfg_file.close()
    return atom_type_table
    ##########################################################################################################################
    
def Generate_SCF(POPULATION,ENS_FOLDER,PREFIX,typical_espresso_header):
    all_scf_files = [os.path.join(ENS_FOLDER + str(POPULATION), f) for f in os.listdir(ENS_FOLDER + str(POPULATION)) if f.startswith("scf_")]
    
    
    # We will generate the input file in a new directory
    
    for file in all_scf_files:
        # Now we are cycling on the scf_ files we found.
        # We must extract the number of the file
        # The file is the string "data_ensemble_manual/scf_population1_X.dat"
        # Therefore the X number is after the last "_" and before the "." character
        # We can split before the string file at each "_", isolate the last part "X.dat"
        # and then split it again on "." (obtaining ["X", "dat"]) and select the first element
        # then we convert the "X" string into an integer
        number = int(file.split("_")[-1].split(".")[0])
        
        # We decide the filename for the espresso input
        # We will call it run_calculation/espresso_run_X.pwi
        filename = os.path.join("scfin", PREFIX + "{}.scf.in".format(number))
        
        # We start writing the file
        with open(filename, "w") as f:
            # We write the header
            f.write(typical_espresso_header)
            
            # Load the scf_population_X.dat file
            ff = open(file, "r")
            structure_lines = ff.readlines()
            ff.close()
            
            # Write the content on the espresso_run_X.pwi file
            # Note in the files we specify the units for both the cell and the structure [Angstrom]
            f.writelines(structure_lines) 
        f.close()   
 
def Generate_POSCAR(POPULATION,ENS_FOLDER,PREFIX):
    all_scf_files = [os.path.join(ENS_FOLDER + str(POPULATION), f) for f in os.listdir(ENS_FOLDER + str(POPULATION)) if f.startswith("scf_")]


    # We will generate the input file in a new directory

    for file in all_scf_files:
        # Now we are cycling on the scf_ files we found.
        # We must extract the number of the file
        # The file is the string "data_ensemble_manual/scf_population1_X.dat"
        # Therefore the X number is after the last "_" and before the "." character
        # We can split before the string file at each "_", isolate the last part "X.dat"
        # and then split it again on "." (obtaining ["X", "dat"]) and select the first element
        # then we convert the "X" string into an integer
        number = int(file.split("_")[-1].split(".")[0])

        # We decide the filename for the espresso input
        # We will call it run_calculation/espresso_run_X.pwi
        filename = os.path.join("scfin", PREFIX + "{}.POSCAR".format(number))

        # We start writing the file
        with open(filename, "w") as f:
            # We write the header

            # Load the scf_population_X.dat file
            ff = open(file, "r")
            structure_lines = ff.readlines()
            ff.close()

            file_lenght = len(structure_lines)
            f.write("\n")
            f.write("1.0000\n")
            for i in range(0,3): f.write(structure_lines[i+1])

            atm_t = []
            atm_n = []
            atom_ordering = [None]*(file_lenght-6)
            atm_t.append(structure_lines[6].split()[0])
            atm_n.append(1)

            for i in range(7,file_lenght):
                flag = True
                for j in range(0,len(atm_t)):
                    if atm_t[j] == structure_lines[i].split()[0]:
                        atm_n[j] = atm_n[j] + 1
                        flag = False
                if flag == True:
                    atm_t.append(structure_lines[i].split()[0])
                    atm_n.append(1)
            for i in atm_t:
                f.write(i + str("  "))
            f.write("\n")
            for i in atm_n:
                f.write(str(i) + str("  "))
            f.write("\nCartesian\n")



            count_n = 0
            for i in atm_t:
                for j in range(6,file_lenght):
                    if i == structure_lines[j].split()[0]:
                        f.write(structure_lines[j].split()[1] + str("  ") + structure_lines[j].split()[2] + str("  ") + structure_lines[j].split()[3] + str("\n"))
                        atom_ordering[j-6] = count_n
                        count_n += 1

            # Write the content on the espresso_run_X.pwi file
            # Note in the files we specify the units for both the cell and the structure [Angstrom]
        f.close()
    return atom_ordering

def read_SCF(POPULATION,ENS_FOLDER,PREFIX,N_RANDOM):
    directory = "scfin"
    output_filenames = [f for f in os.listdir(directory) if f.endswith(".out")] # We select only the output files
    output_files = [os.path.join(directory, f) for f in output_filenames] # We add the directory/outpufilename to load them correctly
    energies = np.zeros(N_RANDOM)
    
    for file in output_files:
        # Get the number of the configuration.
        id_number = int(file.split(".")[0].split(PREFIX)[-1]) # The same as before, we need the to extract the configuration number from the filename
        # Load the file
        ff = open(file, "r")
        lines = [l.strip() for l in ff.readlines()] # Read the whole file removing tailoring spaces
        ff.close()
        
        # Lets look for the energy (in espresso the first line that starts with !)
        # next is used to find only the first occurrence
        energy_line = next(l for l in lines if len(l) > 0 if l.split()[0] == "!")
        
        # Lets collect the energy (the actual number is the 5th item on the line, but python indexes start from 0)
        # note, also the id_number are saved starting from 1
        energies[id_number - 1] = float(energy_line.split()[4])
        
        # Now we can collect the force
        # We need the number of atoms
        nat_line = next( l for l in lines if len(l) > 0 if l.split()[0] == "number" and l.split()[2] == "atoms/cell" )
        nat = int(nat_line.split()[4])
        
        # Now allocate the forces and read them
        forces = np.zeros((nat, 3))
        forces_lines = [l for l in lines if len(l) > 0 if l.split()[0] == "atom"] # All the lines that starts with atom will contain a force
        for i in range(nat):
            forces[i, :] = [float(x) for x in forces_lines[i].split()[-3:]] # Get the last three number from the line containing the force
        
        # Now we can take the stress tensor
        stress = np.zeros((3,3))
        # We pick the index of the line that starts with the words total stress
        index_before_stress = next(i for i, l in enumerate(lines) if len(l) > 0 if l.split()[0] == "total" and l.split()[1] == "stress")
        # The stress tensor is located just after it
        for i in range(3):
            index = i + index_before_stress + 1
            stress[i, :] = [float(x) for x in lines[index].split()[:3]]
    
        # We can save the forces_population1_X.dat and pressures_population1_X.dat files
        force_file = os.path.join( ENS_FOLDER + str(POPULATION), "forces_population"+ str(POPULATION) +"_{}.dat".format(id_number))
        stress_file = os.path.join( ENS_FOLDER + str(POPULATION), "pressures_population" + str(POPULATION) + "_{}.dat".format(id_number))
        np.savetxt(force_file, forces)
        np.savetxt(stress_file, stress)
    
    # Now we read all the configurations, we can save the energy file
    energy_file = os.path.join(ENS_FOLDER + str(POPULATION), "energies_supercell_population" + str(POPULATION) + ".dat")
    np.savetxt(energy_file, energies)
    return energies







def read_OUTCAR(POPULATION,ENS_FOLDER,PREFIX,N_RANDOM,saved_ordering):
    directory = "scfin"
    output_filenames = [f for f in os.listdir(directory) if f.endswith("OUTCAR")] # We select only the output files
    output_files = [os.path.join(directory, f) for f in output_filenames] # We add the directory/outpufilename to load them correctly
    energies = np.zeros(N_RANDOM)

    for file in output_files:
        nat = 0
        # Get the number of the configuration.
        id_number = int(file.split(".")[0].split(PREFIX)[-1]) # The same as before, we need the to extract the configuration number from the filename
        # Load the file
        ff = open(file, "r")
        lines = [l.strip() for l in ff.readlines()] # Read the whole file removing tailoring spaces
        ff.close()

        # Lets look for the energy (in espresso the first line that starts with !)
        # next is used to find only the first occurrence
        energy_line = next(l for l in lines if len(l) > 0 if l.find("free  energy   TOTEN ") != -1)

        # Lets collect the energy (the actual number is the 5th item on the line, but python indexes start from 0)
        # note, also the id_number are saved starting from 1
        energies[id_number - 1] = float(energy_line.split()[4])*0.0734985857

        # Now we can collect the force
        # We need the number of atoms
        nat_line = next( l for l in lines if len(l) > 0 if l.split()[0] == "ions" and l.split()[2] == "type" )


        for i in range(4, len(nat_line.split())):
            nat = nat + int( nat_line.split()[i] )

        # Now allocate the forces and read them
        forces = np.zeros((nat, 3))
        flag = False
        force_lines = []
        for i in range(0,len(lines)):
            if lines[i].find("total drift:") != -1 and flag == True : flag = False
            if flag == True and len(lines[i].split()) > 1: force_lines.append(lines[i])
            if lines[i].find("TOTAL-FORCE (eV/Angst)") != -1 and flag == False : flag = True
            if len(lines[i].split())  > 2:
                if lines[i].split()[0] == "in" and lines[i].split()[1] == "kB": index_before_stress = lines[i]


        print(saved_ordering)
        for i in range(nat):
            forces[i, :] = [float(x)*0.0734985857*0.5029177210544 for x in force_lines[saved_ordering[i]].split()[-3:]] # Get the last three number from the line containing the force


        # Now we can take the stress tensor
        stress = np.zeros((3,3))
        # We pick the index of the line that starts with the words total stress
        # The stress tensor is located just after it
        stress[0,0] = float(index_before_stress.split()[2])/10*math.pow(10,9)/math.pow(10,30)/1.60217663*math.pow(10,19)*0.0734985857*math.pow(0.5291772105,3)
        stress[1,1] = float(index_before_stress.split()[3])/10*math.pow(10,9)/math.pow(10,30)/1.60217663*math.pow(10,19)*0.0734985857*math.pow(0.5291772105,3) 
        stress[2,2] = float(index_before_stress.split()[4])/10*math.pow(10,9)/math.pow(10,30)/1.60217663*math.pow(10,19)*0.0734985857*math.pow(0.5291772105,3)
        stress[0,1] = float(index_before_stress.split()[5])/10*math.pow(10,9)/math.pow(10,30)/1.60217663*math.pow(10,19)*0.0734985857*math.pow(0.5291772105,3)
        stress[1,0] = float(index_before_stress.split()[5])/10*math.pow(10,9)/math.pow(10,30)/1.60217663*math.pow(10,19)*0.0734985857*math.pow(0.5291772105,3)
        stress[0,2] = float(index_before_stress.split()[6])/10*math.pow(10,9)/math.pow(10,30)/1.60217663*math.pow(10,19)*0.0734985857*math.pow(0.5291772105,3)
        stress[2,1] = float(index_before_stress.split()[7])/10*math.pow(10,9)/math.pow(10,30)/1.60217663*math.pow(10,19)*0.0734985857*math.pow(0.5291772105,3)
        stress[1,2] = float(index_before_stress.split()[7])/10*math.pow(10,9)/math.pow(10,30)/1.60217663*math.pow(10,19)*0.0734985857*math.pow(0.5291772105,3)

        # We can save the forces_population1_X.dat and pressures_population1_X.dat files
        force_file = os.path.join( ENS_FOLDER + str(POPULATION), "forces_population"+ str(POPULATION) +"_{}.dat".format(id_number))
        stress_file = os.path.join( ENS_FOLDER + str(POPULATION), "pressures_population" + str(POPULATION) + "_{}.dat".format(id_number))
        np.savetxt(force_file, forces)
        np.savetxt(stress_file, stress)

    # Now we read all the configurations, we can save the energy file
    energy_file = os.path.join(ENS_FOLDER + str(POPULATION), "energies_supercell_population" + str(POPULATION) + ".dat")
    np.savetxt(energy_file, energies)
    return energies















    
    
def Compile_Trainingset(POPULATION,ENS_FOLDER,atom_type_table,conf_table_gamma):
    Ry_to_eV = 13.6057039763
    Bohr_to_A = 0.529177210903

    cfg_file = open("trainset.cfg", "a")
    directory = ENS_FOLDER + str(POPULATION)
    
    energy_file = open(directory + "/energies_supercell_population" + str(POPULATION) + ".dat", "r")
    energy_lines = [l.strip() for l in energy_file.readlines()]
    energy_file.close()
    
    print(len(energy_lines))
    
    for j in range(1,len(conf_table_gamma)+1):
        structure_file = open(directory + "/scf_population" + str(POPULATION) + "_" + str(conf_table_gamma[j-1][0]) + ".dat", "r")
        lines = [l.strip() for l in structure_file.readlines()]
        
        force_file = open(directory + "/forces_population" + str(POPULATION) + "_" + str(conf_table_gamma[j-1][0]) + ".dat", "r")
        force_lines = [l.strip() for l in force_file.readlines()]
        force_file.close()
    
        pressure_file = open(directory + "/pressures_population" + str(POPULATION) + "_" + str(conf_table_gamma[j-1][0]) + ".dat", "r")
        pressure_lines = [l.strip() for l in pressure_file.readlines()]
        pressure_file.close()
        
        SIZE = len(lines)-6
        cell = np.zeros((3, 3))
        atoms = np.zeros((SIZE,3))
        forces = np.zeros((SIZE,3))
        atm_type = [None] * SIZE
        pressures = np.zeros((3, 3))
        
        for i in range(0,3):
            cell[i, :] = [float(x) for x in lines[i+1].split()[-3:]]
            
        for i in range(0,SIZE):
            atoms[i, :] = [float(x) for x in lines[i+6].split()[-3:]]
            atm_type[i] =  lines[i+6].split()[0]
    
            forces[i,:] = [float(x) for x in force_lines[i].split()[-3:]]
            forces[i,:] = forces[i, :] * Ry_to_eV / Bohr_to_A
            
            
        for i in atm_type:
            try: atom_type_table[i]
            except: atom_type_table[i] = len(atom_type_table)    
        
        
        Volume = np.dot(cell[0],np.cross(cell[1], cell[2]))
        
        
        for i in range(0,3):
            pressures[i, :] = [float(x) for x in pressure_lines[i].split()[-3:]]    
            pressures[i, :] = pressures[i, :] * Ry_to_eV / Bohr_to_A / Bohr_to_A / Bohr_to_A* Volume
        
        cfg_file.write("BEGIN_CFG\n")
        cfg_file.write(" Size\n")
        cfg_file.write("{: >5}".format(SIZE) +"\n")
        cfg_file.write(" Supercell\n")
        for row in cell:
            cfg_file.write("    {: >13f} {: >13f} {: >13f}\n".format(*row))
        cfg_file.write(" AtomData:  id type       cartes_x      cartes_y      cartes_z           fx          fy          fz\n")
        for i in range(0,len(atoms)):
            cfg_file.write("    {: >10}".format(i+1) + "{: >5}".format(atom_type_table[atm_type[i]]) + "  {: >13f} {: >13f} {: >13f}".format(*atoms[i,:]) + "  {: >11f} {: >11f} {: >11f}\n".format(*forces[i,:]))
            
        cfg_file.write(" Energy\n")
        cfg_file.write("     {: >13f}".format(float(energy_lines[conf_table_gamma[j-1][0]-1])*Ry_to_eV) +"\n")
        cfg_file.write(" PlusStress:  xx          yy          zz          yz          xz          xy\n")
        cfg_file.write("    {: >12f}".format(pressures[0,0]) + "{: >12f}".format(pressures[1,1]) +  "{: >12f}".format(pressures[2,2]))
        cfg_file.write("{: >12f}".format(pressures[1,2]) + "{: >12f}".format(pressures[0,2]) +  "{: >12f}\n".format(pressures[0,1]))
        cfg_file.write(" Feature atom_type_table " + str(atom_type_table) + "\n")
        cfg_file.write(" Feature conf_number " + str(conf_table_gamma[j-1][0]) + "\n")
        cfg_file.write(" Feature population " + str(POPULATION) + "\n")
        cfg_file.write("END_CFG\n\n")
    cfg_file.close()    
    
def read_CFG(POPULATION,ENS_FOLDER,PREFIX,conf_table_gamma,energies):
    Ry_to_eV = 13.6057039763
    Bohr_to_A = 0.529177210903
    
    directory = ENS_FOLDER + str(POPULATION)
    cfg_file = "MLIP_" + str(POPULATION) + ".out.cfg"
    
    mlip_conf_list = []
    
    for i in range(0,len(conf_table_gamma)):
        mlip_conf_list.append(conf_table_gamma[i][0])
    
    
    conf_output_file = open( cfg_file, "r")
    energy_file = open(directory + "/energies_supercell_population" + str(POPULATION) + ".dat", "w")
    
    ITEM=0
    FORCE_FLAG = False
    ENERGY_FLAG = False
    PRESSURE_FLAG = False
    CELL_FLAG = False
    
    #Volume = np.dot(cell[0],np.cross(cell[1], cell[2]))
    Volume = 1
    
    print(energies)
    print(mlip_conf_list)
    
    TOTAL_NUMBER_OF_MLIP_CONF = 0
    for line in conf_output_file:
    #OPEN FILES ZEROS STORAGE ######################################################################################
        if line.find("BEGIN_CFG") != -1:
            ITEM += 1
            cell = np.zeros((3, 3))
            cell_tmp = []
            forces = []
            pressure = [[0,0,0],[0,0,0],[0,0,0]]
            energy = 0
    #PRINT OUTPUTS CLOSES FILES ######################################################################################
        if line.find("END_CFG") != -1:
            
            
            if conf_number in mlip_conf_list or abs(float(energies[conf_number-1])) > 0.0001:
                
                print("reading " + str(conf_number) + " from DFT" )
                energy_file.write(str(energies[conf_number-1]) + "\n")
                
            else:   
                
                print("reading " + str(ITEM) + " from MLIP" )
                TOTAL_NUMBER_OF_MLIP_CONF = TOTAL_NUMBER_OF_MLIP_CONF +1
                
                pressure_file = open(directory + "/pressures_population" + str(POPULATION) + "_" + str(ITEM) + ".dat", "w")
                force_file = open(directory + "/forces_population" + str(POPULATION) + "_" + str(ITEM) + ".dat", "w")
                
                energy_file.write(str(energy/Ry_to_eV) + "\n")
                
                for data_to_print in forces:
                    force_file.write("{} {} {}".format(data_to_print[0],data_to_print[1],data_to_print[2]) + "\n")
    
                cell[0][:] = cell_tmp[0][:]
                cell[1][:] = cell_tmp[1][:]
                cell[2][:] = cell_tmp[2][:]
                Volume = np.dot(cell[0],np.cross(cell[1], cell[2]))
            
            
                pressure_file.write("{} {} {}".format(pressure[0][0]/Volume,pressure[0][1]/Volume,pressure[0][2]/Volume) + "\n")
                pressure_file.write("{} {} {}".format(pressure[1][0]/Volume,pressure[1][1]/Volume,pressure[1][2]/Volume) + "\n")
                pressure_file.write("{} {} {}".format(pressure[2][0]/Volume,pressure[2][1]/Volume,pressure[2][2]/Volume) + "\n")
          
                pressure_file.close()
                force_file.close()
    #GETS CELL AND STORES ######################################################################################        
        if line.find("Supercell") != -1 or CELL_FLAG == True:
            CELL_FLAG = True
            fields_cell = line.split()
            if fields_cell[0] != "Supercell" and fields_cell[0] != "AtomData:":
                cell_tmp.append([float(fields_cell[0]),float(fields_cell[1]),float(fields_cell[2])])
    #GETS FORCES AND STORES
        if line.find("AtomData:") != -1 or FORCE_FLAG == True:
            CELL_FLAG = False
            FORCE_FLAG = True
            fields_forces=line.split()
            if fields_forces[0] != "AtomData:" and fields_forces[0] != "Energy":
                forces.append([float(fields_forces[5])/ Ry_to_eV * Bohr_to_A,float(fields_forces[6])/ Ry_to_eV * Bohr_to_A,float(fields_forces[7])/ Ry_to_eV * Bohr_to_A])
    #GETS ENERGY AND STORES ######################################################################################
        if line.find("Energy") != -1 or ENERGY_FLAG == True:
            FORCE_FLAG = False
            ENERGY_FLAG = True
            fields_energy=line.split()
            if fields_energy[0] != "Energy" and fields_energy[0] != "PlusStress:":
                energy = float(fields_energy[0])
    #GETS PRESSURES AND STORES ######################################################################################
        if line.find("PlusStress:") != -1 or PRESSURE_FLAG == True:
            ENERGY_FLAG = False
            PRESSURE_FLAG = True
            fields_pressure = line.split()
            if fields_pressure[0] != "PlusStress:" and fields_pressure[0] != "Feature":
                pressure[0][0] = float(fields_pressure[0]) / Ry_to_eV * Bohr_to_A * Bohr_to_A * Bohr_to_A
                pressure[1][1] = float(fields_pressure[1]) / Ry_to_eV * Bohr_to_A * Bohr_to_A * Bohr_to_A
                pressure[2][2] = float(fields_pressure[2]) / Ry_to_eV * Bohr_to_A * Bohr_to_A * Bohr_to_A
                pressure[1][2] = float(fields_pressure[3]) / Ry_to_eV * Bohr_to_A * Bohr_to_A * Bohr_to_A
                pressure[2][1] = float(fields_pressure[3]) / Ry_to_eV * Bohr_to_A * Bohr_to_A * Bohr_to_A
                pressure[0][2] = float(fields_pressure[4]) / Ry_to_eV * Bohr_to_A * Bohr_to_A * Bohr_to_A
                pressure[2][0] = float(fields_pressure[4]) / Ry_to_eV * Bohr_to_A * Bohr_to_A * Bohr_to_A
                pressure[0][1] = float(fields_pressure[5]) / Ry_to_eV * Bohr_to_A * Bohr_to_A * Bohr_to_A
                pressure[1][0] = float(fields_pressure[5]) / Ry_to_eV * Bohr_to_A * Bohr_to_A * Bohr_to_A
    #STOPS LOOP  ######################################################################################
        if line.find("Feature") != -1:
            PRESSURE_FLAG = False
        
        if line.find("Feature   conf_number") != -1:
            conf_number = int(line.split()[-1])  
                                                
    energy_file.close()        
    conf_output_file.close()
    return TOTAL_NUMBER_OF_MLIP_CONF
    
def MTP_CALC_GRADE(POPULATION,MLIP_PATH,PRETRAINED):
#/home/fbelli/Programs/mlip-2-master/bin/mlp
    res = os.system(MLIP_PATH + ' mindist MLIP_' + str(POPULATION) + '.cfg')
    #CALCULATES THE GAMMA FOR THE SELECTION OF THE TRAINING SET
    if POPULATION == 1 and PRETRAINED == False:
        res = os.system(MLIP_PATH + ' calc-grade potential.mtp MLIP_'+str(POPULATION)+'.cfg MLIP_'+str(POPULATION)+'.cfg MLIP_'+str(POPULATION)+'.gamma.cfg')
    else:
        res = os.system(MLIP_PATH + ' calc-grade potential.mtp trainset.cfg MLIP_'+str(POPULATION)+'.cfg MLIP_'+str(POPULATION)+'.gamma.cfg')
 
 
 
def Fill_Gamma_Table(POPULATION,GAMMA,PRETRAINED):
    conf_output_file = open( "MLIP_"+str(POPULATION) + ".gamma.cfg", "r")
    conf_table_gamma = []
             
    for line in conf_output_file:
    #PRINT OUTPUTS CLOSES FILES 
        if line.find("END_CFG") != -1:
            if gamma_conf >= GAMMA or (POPULATION == 1 and PRETRAINED == False):# or abs(float(energy_lines[conf_number-1])) > 0.0001 :
                conf_table_gamma.append([conf_number,gamma_conf])    
    #GETS THE GAMMA FOR ACTIVE LEARNING
        if line.find("Feature   MV_grade") != -1:
            gamma_conf = float(line.split()[-1])
     
        if line.find("Feature   conf_number") != -1:
            conf_number = int(line.split()[-1])   
        
    print(len(conf_table_gamma),conf_table_gamma) 
    conf_output_file.close()
    return conf_table_gamma
    
def Send_to_Folders(PREFIX,LOCAL,PATH,folder,ADDRESS,conf_table_gamma,MAX,marker):
    current_dir = os.getcwd()
    res = os.system("> scfin/list.dat")
    for i in range(1,len(conf_table_gamma)):
        tmp_line = "echo '#" +str(i+1) + "# " + PREFIX + str(conf_table_gamma[i][0]) + "' >> scfin/list.dat"
        res = os.system(tmp_line)   

    Submitter_file = open("Submitter.sh", "w")
    if int(MAX) > int(len(conf_table_gamma)):
        Submitter_file.write('#!/bin/bash -l\n\nmax=' + str(len(conf_table_gamma)) + '\ntop=' + str(len(conf_table_gamma)) +  '\nPREFNAME="' + marker + '" \n')
    else:
        Submitter_file.write('#!/bin/bash -l\n\nmax=' + str(MAX) + '\ntop=' + str(len(conf_table_gamma)) +  '\nPREFNAME="' + marker + '" \n')    
    Submitter_file.close()
    res = os.system('cat base_script.sh >> Submitter.sh' )
    res = os.system('chmod +x Submitter.sh')
        
    if LOCAL == False:    
        res = os.system('scp -q scfin/list.dat Submitter.sh ' + ADDRESS + ':"' + PATH + folder + '/"')
        res = os.system('scp -q scfin/' + PREFIX +'*.in ' + ADDRESS + ':"' + PATH + folder + '/Ins"')
        res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/; ./Submitter.sh"')
    else:
        res = os.system('cp  scfin/list.dat Submitter.sh ' + PATH + folder + '/')
        res = os.system('cp  scfin/' + PREFIX + '*.in ' + PATH + folder + '/Ins')
        res = os.system('cd ' + current_dir)
        print("submitting into send to folder")
        res = os.system('cd ' + PATH + folder + '/; ./Submitter.sh')        
        res = os.system('cd ' + current_dir)
    res = os.system('rm Submitter.sh')


def Send_to_Folders_VASP(PREFIX,LOCAL,PATH,folder,ADDRESS,conf_table_gamma,MAX,marker):
    current_dir = os.getcwd()
    res = os.system("> scfin/list.dat")
    for i in range(0,len(conf_table_gamma)):
        tmp_line = "echo '#" +str(i+1) + "# " + PREFIX + str(conf_table_gamma[i][0]) + "' >> scfin/list.dat"
        res = os.system(tmp_line)

    Submitter_file = open("Submitter.sh", "w")
    if int(MAX) > int(len(conf_table_gamma)):
        Submitter_file.write('#!/bin/bash -l\n\nmax=' + str(len(conf_table_gamma)) + '\ntop=' + str(len(conf_table_gamma)) +  '\nPREFNAME="' + marker + '" \n')
    else:
        Submitter_file.write('#!/bin/bash -l\n\nmax=' + str(MAX) + '\ntop=' + str(len(conf_table_gamma)) +  '\nPREFNAME="' + marker + '" \n')
    Submitter_file.close()
    res = os.system('cat base_script_VASP.sh >> Submitter.sh' )
    res = os.system('chmod +x Submitter.sh')

    if LOCAL == False:
        res = os.system('scp -q scfin/list.dat Submitter.sh ' + ADDRESS + ':"' + PATH + folder + '/"')
        res = os.system('scp -q scfin/' + PREFIX +'*.POSCAR ' + ADDRESS + ':"' + PATH + folder + '/Ins"')
        res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/; ./Submitter.sh"')
    else:
        res = os.system('cp  scfin/list.dat Submitter.sh ' + PATH + folder + '/')
        res = os.system('cp  scfin/' + PREFIX + '*.POSCAR ' + PATH + folder + '/Ins')
        res = os.system('cd ' + current_dir)
        print("submitting into send to folder")
        res = os.system('cd ' + PATH + folder + '/; ./Submitter.sh')
        res = os.system('cd ' + current_dir)
    res = os.system('rm Submitter.sh')




























def Queue_and_resubmit(PREFIX,ADDRESS,PATH,folder,TIMER,LOCAL,conf_table_gamma,MAX,marker):
    current_dir = os.getcwd()
    p=True
   
    while(p==True):
        
        print("Starting cycle.")
        print("Downloading queue info.")
        if LOCAL == False:
            res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/; ./queue_write.sh"')
        else:
            res = os.system('cd ' + current_dir) 
            res = os.system('cd ' + PATH + folder + '/; ./queue_write.sh')
            res = os.system('cd ' + current_dir) 
        print("Downloading output files.")
        if LOCAL == False:

            res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/Outs; tar -czf archive.tar.gz *.out"')
            res = os.system('scp -q ' + ADDRESS + ':"' + PATH + folder + '/Outs/archive.tar.gz" ./scfin/' )
            res = os.system('scp -q ' + ADDRESS + ':"' + PATH + folder + '/queue.tmp" ./')

        else:
            res = os.system('cd ' + current_dir) 
            res = os.system('cd ' + PATH + folder + '/Outs; tar -czf archive.tar.gz *.out')
            res = os.system('cd ' + current_dir) 
            res = os.system('cp ' + PATH + folder + '/Outs/archive.tar.gz ./scfin/' )
            res = os.system('cp ' + PATH + folder + '/queue.tmp ./')
                    
        res = os.system('tar -xf ./scfin/archive.tar.gz -C ./scfin/')       

        queue_file = open("queue.tmp","r")
        queue = queue_file.read()
        print(queue)
        still_in_queue = queue.find(marker)
        
        if still_in_queue != -1:
            print("Still in queue.")
            res = os.system('rm queue.tmp')
            res = os.system('rm ./scfin/archive.tar.gz')
            print("Cleaning on cluster.")

            if LOCAL == False:
                res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/; rm queue.tmp; rm Outs/archive.tar.gz "')
            else:
                res = os.system('cd '+ PATH + folder + '/; rm queue.tmp; rm Outs/archive.tar.gz')           

            ttt = time.localtime()
            print("Standing by from " + str(time.strftime("%H:%M:%S", ttt)) + " for " + str(TIMER/60) + " minutes (" + str(TIMER/3600) + " hours)")
            time.sleep(TIMER)
            
        else:   
            res = os.system('rm queue.tmp')
            res = os.system('rm ./scfin/archive.tar.gz')
            print("Cleaning up on cluster.")
            if LOCAL == False:
                res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/; rm queue.tmp; rm Outs/archive.tar.gz "')
            else:
                res = os.system('cd ' + current_dir)
                res = os.system('cd '+ PATH + folder + '/; rm queue.tmp; rm Outs/archive.tar.gz ')
                res = os.system('cd ' + current_dir)
            
            j=0
            list_file = open("./scfin/list.dat","w")
            for i in range(1,len(conf_table_gamma)+1):
                filename = "./scfin/" + PREFIX + str(conf_table_gamma[i-1][0]) + ".scf.out"
                try:
                    open_output = open(filename, 'r')
                    output_file_stored = open_output.read()
                    check_var = output_file_stored.find("JOB DONE")
                    open_output.close()
                except:
                    check_var = -1
                    
                if check_var == -1:
                    j = j+1
                    list_file.write("#"+str(j)+"# "+ PREFIX + str(conf_table_gamma[i-1][0]) + "\n")
            list_file.close()
            
            if j>0:
                print("Resubmitting jobs.")
                Submitter_in = open("Submitter.sh","w")
                if j < MAX:
                    Submitter_in.write('#!/bin/bash -l\n\nmax=' + str(j) + '\ntop=' + str(j) +  '\nPREFNAME="' + marker + '"\n')
                else:
                    Submitter_in.write('#!/bin/bash -l\n\nmax=' + str(MAX) + '\ntop=' + str(j) +  '\nPREFNAME="' + marker + '"\n')
                Submitter_in.close()
                res = os.system('cat base_script.sh >> Submitter.sh')
                
                if LOCAL == False:    
                    print('scp -q scfin/' + PREFIX + '*.in ' + ADDRESS + ':"' + PATH + folder + '/Ins"')
                    res = os.system('scp -q scfin/list.dat Submitter.sh ' + ADDRESS + ':"' + PATH + folder + '/"')
                    res = os.system('scp -q scfin/' + PREFIX + '*.in ' + ADDRESS + ':"' + PATH + folder + '/Ins"')
                    res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/; ./Submitter.sh"')
                else:
                    res = os.system('cp  scfin/list.dat Submitter.sh ' + PATH + folder + '/')
                    res = os.system('cp  scfin/' + PREFIX +'*.in ' + PATH + folder + '/Ins')
                    res = os.system('cd ' + current_dir) 
                    res = os.system('cd ' + PATH + folder + '/; ./Submitter.sh') 
                    res = os.system('cd ' + current_dir)       
                res = os.system('rm Submitter.sh')
            else:
                p=False
            if j > 0:
                ttt = time.localtime()
                print("Standing by from " + str(time.strftime("%H:%M:%S", ttt)) + " for " + str(TIMER/60) + " minutes (" + str(TIMER/3600) + " hours)")
                time.sleep(TIMER)
    if LOCAL == False:
        res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/; rm Job* *.out; rm Ins/*; rm Outs/*"') 
    else:
        res = os.system('rm '+ PATH + folder + '/Job*; rm '+ PATH + folder + '/*.out; rm '+ PATH + folder + '/Ins/*; rm '+ PATH + folder + '/Outs/*')
        
























def Queue_and_resubmit_VASP(PREFIX,ADDRESS,PATH,folder,TIMER,LOCAL,conf_table_gamma,MAX,marker):
    current_dir = os.getcwd()
    p=True

    while(p==True):

        print("Starting cycle.")
        print("Downloading queue info.")
        if LOCAL == False:
            res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/; ./queue_write.sh"')
        else:
            res = os.system('cd ' + current_dir)
            res = os.system('cd ' + PATH + folder + '/; ./queue_write.sh')
            res = os.system('cd ' + current_dir)
        print("Downloading output files.")
        if LOCAL == False:

            res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/Outs; tar -czf archive.tar.gz *.OUTCAR"')
            res = os.system('scp -q ' + ADDRESS + ':"' + PATH + folder + '/Outs/archive.tar.gz" ./scfin/' )
            res = os.system('scp -q ' + ADDRESS + ':"' + PATH + folder + '/queue.tmp" ./')

        else:
            res = os.system('cd ' + current_dir)
            res = os.system('cd ' + PATH + folder + '/Outs; tar -czf archive.tar.gz *.OUTCAR')
            res = os.system('cd ' + current_dir)
            res = os.system('cp ' + PATH + folder + '/Outs/archive.tar.gz ./scfin/' )
            res = os.system('cp ' + PATH + folder + '/queue.tmp ./')

        res = os.system('tar -xf ./scfin/archive.tar.gz -C ./scfin/')

        queue_file = open("queue.tmp","r")
        queue = queue_file.read()
        print(queue)
        still_in_queue = queue.find(marker)

        if still_in_queue != -1:
            print("Still in queue.")
            res = os.system('rm queue.tmp')
            res = os.system('rm ./scfin/archive.tar.gz')
            print("Cleaning on cluster.")

            if LOCAL == False:
                res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/; rm queue.tmp; rm Outs/archive.tar.gz "')
            else:
                res = os.system('cd '+ PATH + folder + '/; rm queue.tmp; rm Outs/archive.tar.gz')

            ttt = time.localtime()
            print("Standing by from " + str(time.strftime("%H:%M:%S", ttt)) + " for " + str(TIMER/60) + " minutes (" + str(TIMER/3600) + " hours)")
            time.sleep(TIMER)

        else:
            res = os.system('rm queue.tmp')
            res = os.system('rm ./scfin/archive.tar.gz')
            print("Cleaning up on cluster.")
            if LOCAL == False:
                res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/; rm queue.tmp; rm Outs/archive.tar.gz "')
            else:
                res = os.system('cd ' + current_dir)
                res = os.system('cd '+ PATH + folder + '/; rm queue.tmp; rm Outs/archive.tar.gz ')
                res = os.system('cd ' + current_dir)

            j=0
            list_file = open("./scfin/list.dat","w")
            for i in range(1,len(conf_table_gamma)+1):
                filename = "./scfin/" + PREFIX + str(conf_table_gamma[i-1][0]) + ".OUTCAR"
                try:
                    open_output = open(filename, 'r')
                    output_file_stored = open_output.read()
                    check_var = output_file_stored.find("General timing and accounting informations for this job:")
                    open_output.close()
                except:
                    check_var = -1

                if check_var == -1:
                    j = j+1
                    list_file.write("#"+str(j)+"# "+ PREFIX + str(conf_table_gamma[i-1][0]) + "\n")
            list_file.close()

            if j>0:
                print("Resubmitting jobs.")
                Submitter_in = open("Submitter.sh","w")
                if j < MAX:
                    Submitter_in.write('#!/bin/bash -l\n\nmax=' + str(j) + '\ntop=' + str(j) +  '\nPREFNAME="' + marker + '"\n')
                else:
                    Submitter_in.write('#!/bin/bash -l\n\nmax=' + str(MAX) + '\ntop=' + str(j) +  '\nPREFNAME="' + marker + '"\n')
                Submitter_in.close()
                res = os.system('cat base_script_VASP.sh >> Submitter.sh')

                if LOCAL == False:
                    print('scp -q scfin/' + PREFIX + '*.in ' + ADDRESS + ':"' + PATH + folder + '/Ins"')
                    res = os.system('scp -q scfin/list.dat Submitter.sh ' + ADDRESS + ':"' + PATH + folder + '/"')
                    res = os.system('scp -q scfin/' + PREFIX + '*.POSCAR ' + ADDRESS + ':"' + PATH + folder + '/Ins"')
                    res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/; ./Submitter.sh"')
                else:
                    res = os.system('cp  scfin/list.dat Submitter.sh ' + PATH + folder + '/')
                    res = os.system('cp  scfin/' + PREFIX +'*.POSCAR ' + PATH + folder + '/Ins')
                    res = os.system('cd ' + current_dir)
                    res = os.system('cd ' + PATH + folder + '/; ./Submitter.sh')
                    res = os.system('cd ' + current_dir)
                res = os.system('rm Submitter.sh')
            else:
                p=False
            if j > 0:
                ttt = time.localtime()
                print("Standing by from " + str(time.strftime("%H:%M:%S", ttt)) + " for " + str(TIMER/60) + " minutes (" + str(TIMER/3600) + " hours)")
                time.sleep(TIMER)
    if LOCAL == False:
        res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/; rm Job* *.out; rm Ins/*; rm Outs/*"')
    else:
        res = os.system('rm '+ PATH + folder + '/Job*; rm '+ PATH + folder + '/*.out; rm '+ PATH + folder + '/Ins/*; rm '+ PATH + folder + '/Outs/*')
































def Train_submit(POPULATION,ADDRESS,PATH,folder,MLIP_PATH,LOCAL,conf_table_gamma):
    current_dir = os.getcwd()
    if LOCAL == False:
        res = os.system('scp -q potential.mtp trainset.cfg MLIP_' + str(POPULATION) + '.cfg ' + ADDRESS + ':"' + PATH +  folder + '/"' )
    else:
        res = os.system('cp  potential.mtp trainset.cfg MLIP_' + str(POPULATION) + '.cfg ' +  PATH +  folder + '/' )
        
    if len(conf_table_gamma) >= 1:
        print("training")
        if LOCAL == False:
            res = os.system('ssh -q ' + ADDRESS + ' "cd ' + PATH + folder + '/; ./SEND_TRAINING.sh"')
        else:
            res = os.system('cd ' + current_dir) 
            res = os.system('cd ' + PATH + folder + '/; ./SEND_TRAINING.sh') 
            res = os.system('cd ' + current_dir)  
        p=True  
        while(p==True):
            if LOCAL == False:
                res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/;  ./queue_write.sh"')
                res = os.system('scp -q ' + ADDRESS + ':"' + PATH + folder + '/queue.tmp" ./')
                res = os.system('ssh -q ' + ADDRESS + ' "cd '+ PATH + folder + '/; rm queue.tmp"')
            else:
                res = os.system('cd ' + current_dir) 
                res = os.system('cd ' + PATH + folder + '/;  ./queue_write.sh')
                res = os.system('cd ' + current_dir) 
                res = os.system('cp ' + PATH + folder + '/queue.tmp ./')
                res = os.system('rm ' + PATH + folder + '/queue.tmp')            
            queue_file = open("queue.tmp","r")
            queue = queue_file.read()
            still_in_queue = queue.find("test_ta")
            queue_file.close()
            if still_in_queue != -1:
                print("Still in queue.")
                res = os.system('rm queue.tmp')
                ttt = time.localtime()
                print("Standing by from " + str(time.strftime("%H:%M:%S", ttt)) + " for " + str(600/60) + " minutes (" + str(100/3600) + " hours)")
                time.sleep(100)                
            else:
                p=False
    if LOCAL == False:
        res = os.system('scp -q ' + ADDRESS + ':"'+ PATH + folder + '/potential.mtp" ./')
    else:
        res = os.system('cd ' + current_dir)
        res = os.system('cp ' + PATH + folder + '/potential.mtp ./')
    res = os.system(MLIP_PATH + ' calc-efs potential.mtp MLIP_'+str(POPULATION)+'.cfg MLIP_'+str(POPULATION)+'.out.cfg')        
        
