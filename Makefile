NGmerge: NGmerge.c NGmerge.h
	gcc -g -Wall -std=gnu99 -fopenmp -O2 -o NGmerge NGmerge.c -lz
