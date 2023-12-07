#!/bin/bash -l
#
#SBATCH -J <PJ6>
#SBATCH --account=edu23.aqti
#SBATCH -t 10:00
#SBATCH --nodes=1
#SBATCH -p shared
#SBATCH -n 1

#ml gnu7/7.2.0
ml PDC
#ml likwid


#Run your executable
gfortran -o XY.x WXY_OG.f90

srun -n -1 ./XY.x
#for value in 0.5 0.89 2.0; do
 # echo $value > in_${value}.txt
 # echo "in_${value}.txt" | srun -n 1 ./XY.x
#done