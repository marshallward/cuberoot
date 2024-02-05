CC = gcc
FC = gfortran
LD = gfortran
CFLAGS = -g -O3 -mavx -mfma
#FCFLAGS = -fdefault-real-8 -g -O2
FCFLAGS = -fdefault-real-8 -g -O3 -mavx -mfma
LDFLAGS = -lm

#CC = icc
#FC = ifort
#CFLAGS = -g -O2
#CFLAGS = -g -O3 -mavx -mfma
#FCFLAGS = -g -O2 -r8 -fprotect-parens
#FCFLAGS = -g -O2 -r8 -fp-model precise
#FCFLAGS = -g -O3 -mavx -mfma -r8

all: compare timing values test

compare: compare.o cubes.o cbrt_ac.o
	$(LD) $(LDFLAGS) -o $@ $^

timing: timing.o cubes.o cbrt_ac.o
	$(LD) $(LDFLAGS) -o $@ $^

values: values.o cubes.o cbrt_ac.o
	$(LD) $(LDFLAGS) -o $@ $^

test: test.o cubes.o cbrt_ac.o
	$(LD) $(LDFLAGS) -o $@ $^

compare.o: compare.f90 cubes.o
	$(FC) $(FCFLAGS) -c -o $@ $<

values.o: values.f90 cubes.o
	$(FC) $(FCFLAGS) -c -o $@ $<

timing.o: timing.f90 cubes.o
	$(FC) $(FCFLAGS) -c -o $@ $<

test.o: test.f90 cubes.o
	$(FC) $(FCFLAGS) -c -o $@ $<

cubes.o: cubes.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

cbrt_ac.o: cbrt_ac.c
	$(CC) $(CFLAGS) -c -o $@ $<

# Plot generation

err.svg: plot.py err_quad.txt
	python plot.py

err_quad.txt: compare
	./compare

clean:
	$(RM) cubes.mod cubes.o timing compare values test
