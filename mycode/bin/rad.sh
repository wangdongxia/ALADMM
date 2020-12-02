#!/bin/sh
#BSUB -e /data/home/xjsjleiyongmei/wdx/RADADMM/paper.err
#BSUB -o /data/home/xjsjleiyongmei/wdx/RADADMM/paper65_4.out
#BSUB -n 65
#BSUB -q priority
#BSUB -J admm
#BSUB -R "span[ptile=8]"

#export LD_LIBRARY_PATH=/data/home/xjsjleiyongmei/xjy/mpich/lib:$LD_LIBRARY_PATH

HOME=/data/home/xjsjleiyongmei/wdx

 



mpirun -machine $LSB_DJOB_HOSTFILE -np 65 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 0  64
mpirun -machine $LSB_DJOB_HOSTFILE -np 65 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 1  64
mpirun -machine $LSB_DJOB_HOSTFILE -np 65 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 2  64
mpirun -machine $LSB_DJOB_HOSTFILE -np 65 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 5  64  
 
mpirun -machine $LSB_DJOB_HOSTFILE -np 65 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 0  16
mpirun -machine $LSB_DJOB_HOSTFILE -np 65 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 1  16
  
 
mpirun -machine $LSB_DJOB_HOSTFILE -np 65 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 0  32
mpirun -machine $LSB_DJOB_HOSTFILE -np 65 $HOME/RADADMM/bin/admm -file $HOME/RADADMM/admm.conf 1  32
  


  
