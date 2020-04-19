#!/bin/sh
#QSUB -queue i18cpu
#QSUB -node 4
#QSUB -mpi 96
#QSUB -omp 1
#QSUB -place pack
#QSUB -over false

#PBS -l walltime=0:30:00
#PBS -N Ar-GK

cd $PBS_O_WORKDIR
nprocs=120
. /etc/profile.d/modules.sh
module load intel/16.0.4.258
module load intel-mpi/5.1.3.258
#module load pgi
EXE=~/lmp_mpi

for i in {1..10} ; do

  dir_name="seed${i}"
  mkdir $dir_name
  v_seed="$i$i$i"
  cp Ar_structure.lmp in.GK  $dir_name
  cd $dir_name
  pwd
  sed -i "s/variable seed.*/variable seed equal $v_seed/g" in.GK
  mpijob   $EXE < in.GK > out
  cd ..
done
wait
