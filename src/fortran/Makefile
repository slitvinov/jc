M=ode
O=ode.o fcn.o
FC=gfortran
FCFLAGS=-O2 -Wall -Wextra

$M: $O;  $(FC) -o $@ $^

ode.o: ode.f90; $(FC) $(FCFLAGS) -c $<
fcn.o: fcn.f90; $(FC) $(FCFLAGS) -c $<

test: $M; atest test/*

clean:; rm -f $M $O

.PHONY: clean test
.SUFFIXES:
