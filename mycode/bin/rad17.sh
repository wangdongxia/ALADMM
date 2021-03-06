#!/bin/sh
#BSUB -e /data/home/xjsjleiyongmei/wdx/RADADMM/paper.err
#BSUB -o /data/home/xjsjleiyongmei/wdx/RADADMM/paper17.out
#BSUB -n 17
#BSUB -q priority
#BSUB -J admm
#BSUB -R "span[ptile=4]"

#export LD_LIBRARY_PATH=/data/home/xjsjleiyongmei/xjy/mpich/lib:$LD_LIBRARY_PATH

HOME=/data/home/xjsjleiyongmei/wdx

 
 
 
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 0  16
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 1  16
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 2  16
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 3  16
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 4  16
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 5  16
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 6  16
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 0  8
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 1  8
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 2  8
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 3  8
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 4  8
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 5  8
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 6  8
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 0  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 1  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 2  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 3  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 4  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 5  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 6  4
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 0  2
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 1  2
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 2  2
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 3  2
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 4  2
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 5  2
mpirun -machine $LSB_DJOB_HOSTFILE -np 17 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 6  2

 


  
