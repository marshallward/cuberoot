FC = gfortran
FCFLAGS = -g -O3 -mavx -mfma -fdefault-real-8
#FCFLAGS = -g -O2 -fdefault-real-8

all: compare timing

compare: compare.f90
	$(FC) $(FCFLAGS) -o $@ $<

timing: timing.f90
	$(FC) $(FCFLAGS) -o $@ $<
