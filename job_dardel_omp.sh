#!/bin/bash -l
#
#SBATCH -J <PJ6>
#SBATCH --account=edu23.aqti
#SBATCH -t 10:00
#SBATCH --nodes=1
#SBATCH -p shared
#SBATCH -n 10


#ml gnu7/7.2.0
ml PDC
#ml likwid


#Run your executable
gfortran -o XY_OMP.x -fopenmp WXY_OMP.f90
#gfortran -o XY_SIMD.x WXY_simd.f90

#num_iter=$(awk '/INTEGER, PARAMETER :: Niter/ {print $NF}' WXY_OMP.f90)
#export OMP_NUM_THREADS=1
#sed -i "s/num/$num_iter/" job_DARDEL_OMP.sh
./XY_OMP.x
#srun -n 1 -c 4 ./XY_SIMD.x
#for value in 0.5 0.89 2.0; do
 # echo $value > in_${value}.txt
 # echo "in_${value}.txt" | srun -n 1 ./XY.x
#done
