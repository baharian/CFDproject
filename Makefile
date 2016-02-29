# Makefile for program #5             ATMS 502/CSE 566, Fall 2012
#
# Lines starting with "#" are comments.
# First line is default, so typing "make" makes executable named pgm4.
# Typing "make test_interp" compiles and creates executable test_interp.
#
# The executable has its name, then all dependencies (object files).
#   Beneath that is the statement to link them and create the program.
#
# The last statements say how to turn .f90 or .c files into .o files: compiling.
#

CC = icc
OPTIONS = -O2 -openmp -fno-alias -fno-fnalias -fargument-noalias -restrict
OBJECTS = putfield.o ic.o bc.o advect.o advect1D.o diffuse.o pgf.o update.o stats.o p5.o

p5:	$(OBJECTS)
	$(CC) $(OPTIONS) -o p5 $(OBJECTS)

.c.o:
	$(CC) $(OPTIONS) -c $*.c
