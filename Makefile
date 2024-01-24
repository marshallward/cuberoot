#FC = gfortran
#FCFLAGS = -g -O2 -fdefault-real-8
#FCFLAGS = -g -O3 -mavx -mfma -fdefault-real-8

FC = ifort
FCFLAGS = -g -O3 -mavx -mfma -r8
#FCFLAGS = -g -O2 -r8

all: compare timing

compare: compare.f90 cubes.o
	$(FC) $(FCFLAGS) -o $@ $^

timing: timing.f90 cubes.o
	$(FC) $(FCFLAGS) -o $@ $^

cubes.o: cubes.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

clean:
	$(RM) cubes.mod cubes.o timing compare
