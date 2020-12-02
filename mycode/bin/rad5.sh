#!/bin/sh
#BSUB -e /data/home/xjsjleiyongmei/wdx/RADADMM/paper.err
#BSUB -o /data/home/xjsjleiyongmei/wdx/RADADMM/paper5.out
#BSUB -n 5
#BSUB -q priority
#BSUB -J admm
#BSUB -R "span[ptile=1]"

#export LD_LIBRARY_PATH=/data/home/xjsjleiyongmei/xjy/mpich/lib:$LD_LIBRARY_PATH

HOME=/data/home/xjsjleiyongmei/wdx


 
 
 


mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 0  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 1  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 2  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 3  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 4  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 5  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 6  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 0  2
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 1  2
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 2  2
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 3  2
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 4  2
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 5  2
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 6  2
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 0  1
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 1  1
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 2  1
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 3  1
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 4  1
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 5  1
mpirun -machine $LSB_DJOB_HOSTFILE -np 5 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 6  1
 

 


  
