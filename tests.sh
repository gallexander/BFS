#!/bin/bash

mpicc -Wall -o gen kronecker_generator_par.c -lm
mpirun -np 4 gen

