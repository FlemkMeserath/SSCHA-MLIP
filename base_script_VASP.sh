
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
#SBATCH --account=ezurek
#SBATCH --clusters=faculty --partition=scavenger --qos=scavenger
##SBATCH --mail-user=someone@buffalo.edu --mail-type=END
##SBATCH -F /projects/academic/ezurek/masashik/templates/nodefile.txt

##########################################
echo '---- Loading Dependencies ----'
module purge
module load intel
ulimit -s unlimited
export I_MPI_PMI_LIBRARY=/opt/software/slurm/lib64/libpmi.so
export vasp="/projects/academic/ezurek/software/vasp6.4.2/vasp.6.4.2/bin/vasp_std"

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

cp -r \$SLURM_SUBMIT_DIR/POTCAR \$WORKDIR
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

