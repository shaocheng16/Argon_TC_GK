#!/bin/sh
#QSUB -queue F36cpu
#QSUB -node  10
#QSUB -place pack
#QSUB -over false

#PBS -l walltime=20:30:00
#PBS -N PbS-GK

cd $PBS_O_WORKDIR
. /etc/profile.d/modules.sh
module load intel/16.0.4.258
module load intel-mpi/5.1.3.258
#module load pgi
EXE=~/lmp_mpi

lmp_input=in.GK
lmp_dat=data.PbS_ligand_SL.lmp
pot_file=in.potential

for i in {1..10} ; do

  dir_name="seed${i}"
  mkdir $dir_name
  v_seed="$i$i$i"
  cp  ${lmp_input} ${lmp_dat} ${pot_file} $dir_name
  cd $dir_name
  pwd
  sed -i "s/variable seed.*/variable seed equal $v_seed/g" ${lmp_input}
  cd ..
done
bulkjob ./job_list
