#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=24 --mem=0
#SBATCH --time=48:00:00
#SBATCH --job-name=test_tail --output=train.res
#SBATCH --account=ezurek
#SBATCH --clusters=faculty --partition=ezurek --qos=ezurek
#SBATCH --requeue
#SBATCH --exclude=cpn-v11-15-02

echo "Loading Modules"
module load  foss

echo '---- Initial Time ----'
echo 'Current time is:'
date
echo "I got here"

#srun -n 24 /projects/academic/ezurek/francesco/mlip-2-master/bin/mlp calc-grade potential.mtp trainset.cfg MLIP_1.cfg MLIP_1.gamma.cfg

mpirun -n 24 /projects/academic/ezurek/francesco/mlip-2-master/bin/mlp train potential.mtp trainset.cfg  --max-iter=600  --trained-pot-name=potential.mtp

#mpirun -n 24 mlp calc-efs potential.mtp MLIP_5.cfg MLIP_5.out.cfg

# mlp calc-errors pot_${lev}_${cmpd}_frx.mtp validate_frx.cfg >> err_${lev}_frx.txt
# mlp calc-efs pot_${lev}_${cmpd}_frx.mtp validate_frx.cfg post_validate_${lev}.cfg
# python ../../get_efs.py validate_frx.cfg post_validate_${lev}.cfg
# tail post_validate_${lev}.csv >> err_${lev}_frx.txt
# rm post_validate*
#done

#for i in 10 16 22
#do 
# cd $pth/$cmpd
# cp ../${i}.mtp ./.
# mind=$(mlp mindist aflow_MTP_${cmpd}.cfg)
# python ../insert_mindist.py ${i}.mtp $mind
# python ../test_train.py aflow_MTP_${cmpd}.cfg train_frx.cfg test_frx.cfg
# srun --propagate=STACK mlp train ${i}.mtp train_frx.cfg --trained-pot-name=pot_${i}_${cmpd}_frx.mtp --valid-cfgs=test_frx.cfg
# rm ${i}.mtp
#done

echo '---- MLIP Job Done ----'
echo 'Current time is:'
date
           
