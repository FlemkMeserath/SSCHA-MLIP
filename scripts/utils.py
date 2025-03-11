import numpy as np
import sys,os
import time
import spglib as spg
import math
from ase.build import bulk


INCAR_FILE = """# output options
LWAVE  = .FALSE. # write or don't write WAVECAR
LCHARG = .FALSE. # write or don't write CHG and CHGCAR
LELF   = .FALSE. # write ELF
LASPH = .TRUE.
ISTART = 0
NSW = 0        # number of ionic steps
IBRION = -1       # 2=conjucate gradient, 1=Newton like
EDIFF = 1E-8     # 1E-3 very low precision for pre-relaxation, use 1E-5 next
PREC = accurate  # precision low, med, high, accurate
KSPACING = 0.1047198
ISMEAR = 0       # -5 = tetraedon, 1..N = Methfessel
SIGMA = 0.02
ENCUT = 1000      # cutoff energy
PSTRESS = 390
NWRITE = 0
#METAGGA = R2SCAN
NCORE = 4
ADDGRID = .TRUE.
ISYM = 0
LREAL=.FALSE.
NWRITE = 0
"""

SCF_FILE = """&control
    calculation        = 'scf'
    restart_mode       = 'from_scratch'
    prefix             = 'scf'
    tstress            = .true. 
    tprnfor            = .true.
    pseudo_dir         = './Pseudo'
    outdir             = './tmp/'
    verbosity          = 'high'
&end

&system
    ibrav              = 0
    nat                = 32 
    ntyp               = 3
    ecutwfc            = 90
    ecutrho            = 900
    degauss            = 0.02
    occupations        = 'smearing'
    smearing           = 'mp'
&end

&electrons
   mixing_beta = 0.2
      conv_thr = 1.0d-9
&end

&ions
&end

ATOMIC_SPECIES

    Pd 106.42   Pd.pbe-n-rrkjus_psl.1.0.0.UPF
    Cu 63.546   Cu.pbe-dn-rrkjus_psl.1.0.0.UPF 
    H  1.007    H.pbe-rrkjus_psl.1.0.0.UPF
 
K_POINTS automatic

 10 10 7 1 1 1 
"""


base_script_QE = """
NODES=1
PROCS=24
POOLS=1
CTIME="24:00:00"

LISTFILE="list.dat"
COMNAME="pw.x"

for i in `seq 1 1 $max`
do

j=$((i+skip))

echo "#!/bin/bash

#SBATCH --nodes=${NODES} --ntasks-per-node=${PROCS} --mem=0
#SBATCH --time=${CTIME}
#SBATCH --job-name=${PREFNAME}${i} --output=${PREFNAME}${i}.out
[ADD HERE THE RESTO OF YOUR SBATCH DETAILS]


##########################################
echo '---- Loading Dependencies ----'
module purge
[load your modules here]


#############################################

" > Job$i.sh

for k in `seq $i $max $top`
do



 SW=`cat $LISTFILE | grep "#$k#"`
 NB=`echo $SW | cut -d " " -f 2`

if [ -z "$NB" ]
then
echo "Cannot find index $NB"
else
echo "
#--------------------------------------------------------------

mkdir \$SLURM_SUBMIT_DIR/$NB
WORKDIR=\${SLURM_SUBMIT_DIR}/$NB

cp -r \$SLURM_SUBMIT_DIR/Pseudo \$WORKDIR
cp \$SLURM_SUBMIT_DIR/Ins/$NB.scf.in \$WORKDIR

cd \$WORKDIR

srun -n \$SLURM_NPROCS [path to executable]/pw.x < $NB.scf.in > $NB.scf.out

cp \$WORKDIR/$NB.scf.out \$SLURM_SUBMIT_DIR/Outs

cd \$SLURM_SUBMIT_DIR/

rm -r \$WORKDIR


#-------------------------------------------------------------
" >> Job$i.sh
fi

done
done


for i in `seq 1 1 $max`
do
sbatch Job$i.sh
done"""


base_script_VASP = """
NODES=1
PROCS=24
POOLS=1
CTIME="24:00:00"

LISTFILE="list.dat"
COMNAME="pw.x"

for i in `seq 1 1 $max`
do

j=$((i+skip))

echo "#!/bin/bash

#SBATCH --nodes=${NODES} --ntasks-per-node=${PROCS} --mem=0
#SBATCH --time=${CTIME}
#SBATCH --job-name=${PREFNAME}${i} --output=${PREFNAME}${i}.out
[ADD HERE THE RESTO OF YOUR SBATCH DETAILS]


##########################################
echo '---- Loading Dependencies ----'
module purge
[load your modules here]
vasp=[place your path to vasp here]


#############################################

" > Job$i.sh

for k in `seq $i $max $top`
do



 SW=`cat $LISTFILE | grep "#$k#"`
 NB=`echo $SW | cut -d " " -f 2`

if [ -z "$NB" ]
then
echo "Cannot find index $NB"
else
echo "
#--------------------------------------------------------------

mkdir \$SLURM_SUBMIT_DIR/$NB
WORKDIR=\${SLURM_SUBMIT_DIR}/$NB

cp -r \$SLURM_SUBMIT_DIR/Pseudo/POTCAR \$WORKDIR
cp -r \$SLURM_SUBMIT_DIR/INCAR \$WORKDIR
cp -r \$SLURM_SUBMIT_DIR/KPOINTS \$WORKDIR
cp \$SLURM_SUBMIT_DIR/Ins/$NB.POSCAR \$WORKDIR/POSCAR

cd \$WORKDIR

srun \$vasp 

cp \$WORKDIR/OUTCAR \$SLURM_SUBMIT_DIR/Outs/$NB.OUTCAR

cd \$SLURM_SUBMIT_DIR/

rm -r \$WORKDIR


#-------------------------------------------------------------
" >> Job$i.sh
fi

done
done


for i in `seq 1 1 $max`
do
sbatch Job$i.sh
done
"""


queue_write = """#!/bin/bash -l
squeue -u fbelli -M faculty > queue.tmp"""

send_training = """#!/bin/bash -l
sbatch train.sh"""


training = """#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=24 --mem=0
#SBATCH --time=48:00:00
#SBATCH --job-name=test_tail --output=train.res
[place here the rest of your sbatch details]

echo "Loading Modules"
[load the rest of your modules here]

echo '---- Initial Time ----'
echo 'Current time is:'
date


mpirun -n 24 [path to executable]/mlp train potential.mtp trainset.cfg  --max-iter=600  --trained-pot-name=potential.mtp

echo '---- MLIP Job Done ----'
echo 'Current time is:'
date"""




def create_scripts(path,CALCULATOR):
    if not os.path.exists( path):
        os.mkdir(path)
    if not os.path.exists(path + "/Ins"):
        os.mkdir( path + "/Ins")
    if not os.path.exists(path + "/Outs"):
        os.mkdir( path + "/Outs")
    res = os.system("cp -r Pseudo " + path + "/")

    file_script = open(path + "/queue_write.sh","w")
    file_script.write(queue_write)
    file_script.close()
    res = os.system("chmod +x "+ path + "/queue_write.sh")

    file_script = open(path + "/SEND_TRAINING.sh","w")
    file_script.write(send_training)
    file_script.close()
    res = os.system("chmod +x "+ path + "/SEND_TRAINING.sh")

    file_script = open(path + "/train.sh","w")
    file_script.write(training)
    file_script.close()
    res = os.system("chmod +x "+ path + "/train.sh")

    if CALCULATOR == "VASP":
        file_script = open(path + "/INCAR","w")
        file_script.write(INCAR_FILE)
        file_script.close()

    
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
    
    header = """&control
        calculation        = 'scf'
        restart_mode       = 'from_scratch'
        prefix             = 'scf'
        tstress            = .true. 
        tprnfor            = .true.
        pseudo_dir         = './Pseudo'
        outdir             = './tmp/'
        verbosity          = 'high'
    &end
    
    &system
        ibrav              = 0
        nat                = 216
        ntyp               = 1
        ecutwfc            = 70.0
        ecutrho            = 700
        degauss            = 0.01
        occupations        = 'smearing'
        smearing           = 'mp'
     &end
    
    &electrons
        mixing_beta        = 0.2
        conv_thr           = 1.0d-9
     &end
    &ions
     &end
    
     ATOMIC_SPECIES
      Li  6.941   Li.pbe-sl-rrkjus_psl.1.0.0.UPF
    
    K_POINTS (automatic)
      3 3 3   1  1  1
"""

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







