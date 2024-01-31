FC = gfortran
#FCFLAGS = -fdefault-real-8 -g -O2
FCFLAGS = -fdefault-real-8 -g -O3 -mavx -mfma

#FC = ifort
#FCFLAGS = -g -O2 -r8 -fprotect-parens
#FCFLAGS = -g -O2 -r8 -fp-model precise
#FCFLAGS = -g -O3 -mavx -mfma -r8

all: compare timing values

compare: compare.f90 cubes.o
	$(FC) $(FCFLAGS) -o $@ $^

values: values.f90 cubes.o
	$(FC) $(FCFLAGS) -o $@ $^

timing: timing.f90 cubes.o
	$(FC) $(FCFLAGS) -o $@ $^

test: test.f90 cubes.o
	$(FC) $(FCFLAGS) -o $@ $^

cubes.o: cubes.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

# Plot generation

err.svg: plot.py err_quad.txt
	python plot.py

err_quad.txt: compare
	./compare

clean:
	$(RM) cubes.mod cubes.o timing compare values test
