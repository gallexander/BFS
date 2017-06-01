#!/bin/bash

PATH=$PATH":"/opt/slurm/bin/
module load mpi/openmpiS

mpicc -Wall kronecker_generator_par.c -o gen -lm
srun -p q_student -N 32 --ntasks-per-node=8 ./gen

mpicc -Wall bfs_seq.c -o bfs_seq -lm
srun -p q_student -N 32 --ntasks-per-node=8 ./bfs_seq

#mpicc -Wall bfs_par_allvisited_parallelsort_improvement.c -o bfs_allreduce_normal -lm
#srun -p q_student -N 32 --ntasks-per-node=1 ./bfs_allreduce_normal
#srun -p q_student -N 32 --ntasks-per-node=2 ./bfs_allreduce_normal
#srun -p q_student -N 32 --ntasks-per-node=4 ./bfs_allreduce_normal
#srun -p q_student -N 32 --ntasks-per-node=8 ./bfs_allreduce_normal
#srun -p q_student -N 32 --ntasks-per-node=16 ./bfs_allreduce_normal

#mpicc -Wall -fopenmp bfs_hybrid_improvement.c -o bfs_hybrid -lm
#export OMP_NUM_THREADS=1
#srun -p q_student -N 32 --ntasks-per-node=1 -c 1 ./bfs_hybrid
#export OMP_NUM_THREADS=2
#srun -p q_student -N 32 --ntasks-per-node=1 -c 2 ./bfs_hybrid
#export OMP_NUM_THREADS=4
#srun -p q_student -N 32 --ntasks-per-node=1 -c 4 ./bfs_hybrid
#export OMP_NUM_THREADS=8
#srun -p q_student -N 32 --ntasks-per-node=1 -c 8 ./bfs_hybrid
#export OMP_NUM_THREADS=16
#srun -p q_student -N 32 --ntasks-per-node=1 -c 16 ./bfs_hybrid

#export OMP_NUM_THREADS=4
#srun -p q_student -N 32 --ntasks-per-node=2 -c 4 ./bfs_hybrid
#export OMP_NUM_THREADS=4
#srun -p q_student -N 32 --ntasks-per-node=4 -c 4 ./bfs_hybrid

mpicc -Wall bfs_par_alltoall.c -o bfs_alltoall -lm
srun -p q_student -N 32 --ntasks-per-node=1 ./bfs_alltoall
srun -p q_student -N 32 --ntasks-per-node=2 ./bfs_alltoall
srun -p q_student -N 32 --ntasks-per-node=4 ./bfs_alltoall
srun -p q_student -N 32 --ntasks-per-node=8 ./bfs_alltoall
srun -p q_student -N 32 --ntasks-per-node=16 ./bfs_alltoall

#mpicc -Wall bfs_par_reducescatter.c -o bfs_reducescatter -lm
#srun -p q_student -N 32 --ntasks-per-node=1 ./bfs_reducescatter
#srun -p q_student -N 32 --ntasks-per-node=2 ./bfs_reducescatter
#srun -p q_student -N 32 --ntasks-per-node=4 ./bfs_reducescatter
#srun -p q_student -N 32 --ntasks-per-node=8 ./bfs_reducescatter
#srun -p q_student -N 32 --ntasks-per-node=16 ./bfs_reducescatter


