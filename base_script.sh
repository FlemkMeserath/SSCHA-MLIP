
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
#SBATCH -F /projects/academic/ezurek/xiaoyu/od/nodefile.txt


##########################################
echo '---- Loading Dependencies ----'
module purge
ulimit -s unlimited
export I_MPI_PMI_LIBRARY=/opt/software/slurm/lib64/libpmi.so
module load intel


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

srun -n \$SLURM_NPROCS /projects/academic/ezurek/software/qe-7.3-sscha-intel/bin/pw.x < $NB.scf.in > $NB.scf.out

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
done

