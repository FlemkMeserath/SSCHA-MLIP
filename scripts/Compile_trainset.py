import sys
import numpy as np

def Compile_Trainingset(trainset_file,ENS_FOLDER,POPULATION):
    __Ry_to_eV__  =  13.6057039763 
    __Bohr_to_A__ = 0.529177249
    atom_type_table = {} 
 
 
    cfg_file = open(trainset_file, "a")
    directory = ENS_FOLDER 
 
    energy_file = open(directory + "/energies_supercell_population" + str(POPULATION) + ".dat", "r")
    energy_lines = [l.strip() for l in energy_file.readlines()]
    energy_file.close()
 
    print(len(energy_lines))
 
    for j in range(1,len(energy_lines)+1):
        structure_file = open(directory + "/scf_population" + str(POPULATION) + "_" + str(j) + ".dat", "r")
        lines = [l.strip() for l in structure_file.readlines()]
 
        force_file = open(directory + "/forces_population" + str(POPULATION) + "_" + str(j) + ".dat", "r")
        force_lines = [l.strip() for l in force_file.readlines()]
        force_file.close()
 
        pressure_file = open(directory + "/pressures_population" + str(POPULATION) + "_" + str(j) + ".dat", "r")
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
            forces[i,:] = forces[i, :] * __Ry_to_eV__ / __Bohr_to_A__
 
 
        for i in atm_type:
            try: atom_type_table[i]
            except: atom_type_table[i] = len(atom_type_table)
 
 
        Volume = np.dot(cell[0],np.cross(cell[1], cell[2]))
 
 
        for i in range(0,3):
            pressures[i, :] = [float(x) for x in pressure_lines[i].split()[-3:]]
            pressures[i, :] = pressures[i, :] * __Ry_to_eV__ / __Bohr_to_A__ / __Bohr_to_A__ / __Bohr_to_A__* Volume
 
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
        cfg_file.write("     {: >13f}".format(float(energy_lines[j-1])*__Ry_to_eV__) +"\n")
        cfg_file.write(" PlusStress:  xx          yy          zz          yz          xz          xy\n")
        cfg_file.write("    {: >12f}".format(pressures[0,0]) + "{: >12f}".format(pressures[1,1]) +  "{: >12f}".format(pressures[2,2]))
        cfg_file.write("{: >12f}".format(pressures[1,2]) + "{: >12f}".format(pressures[0,2]) +  "{: >12f}\n".format(pressures[0,1]))
        cfg_file.write(" Feature atom_type_table " + str(atom_type_table) + "\n")
        cfg_file.write(" Feature conf_number " + str(j) + "\n")
        cfg_file.write(" Feature population " + str(POPULATION) + "\n")
        cfg_file.write("END_CFG\n\n")
    cfg_file.close()


if len(sys.argv) < 4:
    print("\n----------------------------------------------------------------")
    print("This script generates a CFG file for MLIP from the SSCHA EFS files.")
    print("Provide in order: [CFG file location/name] [SSCHA EFS folder] [SSCHA Population]\n")
else:
    input_vars = sys.argv[1], sys.argv[2], sys.argv[3]
    Compile_Trainingset(sys.argv[1], sys.argv[2], int(sys.argv[3]))

