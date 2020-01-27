#
# Makefile for w2_3d_main.c RUN: typing make in this directory
#
# CFLAGS=-c -O3 -Wall -g
CC=gcc
#CC=icc
# CFLAGS= -O3 -Wall -ggdb3
 CFLAGS= -O3 -Wall -std=c99
# $(CC) $(CFLAGS)

all: w2_exec w2_main.o adisco_1d.o inidat.o  integ.o rkc.o w2_eq.o derivs_1d.o binio.o util.o input.o pygraph.o

w2_exec: w2_main.c first_macro_1d.h structs_1d.h gen_1d.h w2_main.o inidat.o integ.o rkc.o adisco_1d.o w2_eq.o derivs_1d.o binio.o util.o input.o pygraph.o
	$(CC) $(CFLAGS) -L/usr/local/lib -L. -I/usr/include -I/usr/include/hdf5/serial -ffast-math -finline-functions -fopenmp -o w2_exec w2_main.o inidat.o integ.o rkc.o adisco_1d.o w2_eq.o derivs_1d.o binio.o util.o input.o pygraph.o  -lm -lhdf5

w2_main.o: first_macro_1d.h structs_1d.h gen_1d.h w2_main.c
	$(CC) $(CFLAGS) -c -o w2_main.o w2_main.c

adisco_1d.o: first_macro_1d.h structs_1d.h gen_1d.h adisco_1d.c
	$(CC) $(CFLAGS) -c -o adisco_1d.o adisco_1d.c

inidat.o: first_macro_1d.h structs_1d.h gen_1d.h  inidat.c
	$(CC) $(CFLAGS) -c -o inidat.o inidat.c

integ.o: first_macro_1d.h structs_1d.h gen_1d.h integ.c
	$(CC) $(CFLAGS) -c -finline-functions -fopenmp -o integ.o integ.c

rkc.o: first_macro_1d.h structs_1d.h gen_1d.h rkc.c
	$(CC) $(CFLAGS) -c -finline-functions -o rkc.o rkc.c

w2_eq.o: w2_eq.c first_macro_1d.h structs_1d.h gen_1d.h
	$(CC) $(CFLAGS) -c -finline-functions -o w2_eq.o w2_eq.c

input.o: input.c first_macro_1d.h structs_1d.h gen_1d.h
	$(CC) $(CFLAGS) -c -o input.o input.c

#output.o: output.c first_macro_1d.h structs_1d.h gen_1d.h
#	$(CC) $(CFLAGS) -c -o output.o output.c

derivs_1d.o: derivs_1d.c first_macro_1d.h structs_1d.h gen_1d.h derivs_1d.h
	$(CC) $(CFLAGS) -c -ffast-math -finline-functions -o derivs_1d.o derivs_1d.c

binio.o: binio.c
	$(CC) $(CFLAGS) -c -DLITTLE -o binio.o binio.c

util.o: util.c
	$(CC) $(CFLAGS) -c -o util.o util.c

pygraph.o: pygraph.c pygraph.h
	$(CC) $(CFLAGS) -I/usr/lib/openmpi/include -I/usr/include -I/usr/include/hdf5/serial -c -o pygraph.o pygraph.c -lhdf5

clean:
	rm -f *.o w2_1d_exec tranSVFile core *~

limpia: 
	rm -f Dump/run_*
