# -*- Makefile -*-

default: quadpack.o quadrule.o rfsfns.o

%.o: %.f90
	$(FC) $(FFLAGS) $(PIXFLAGS) -c -o $@ $<

%.o: %.f
	$(FC) $(FFLAGS) $(PIXFLAGS) -c -o $@ $<

%.o: %.F90
	@echo compiling $*.F90
	@case $(OS) in \
	Darwin ) \
	ln -s -f $< $*.junk.F90; \
	$(FC) $(FFLAGS) $(ADDUS) $(PPFLAGS) -c $*.junk.F90 -o $@; \
	rm $*.junk.F90;; \
	OSF ) \
	ln -s -f $< $*.f90; \
	$(FC) $(FFLAGS) $(ADDUS) -cpp -c $(basename $<).f90 -o $@;; \
	* ) \
	$(FC) $(FFLAGS) $(ADDUS) $(PPFLAGS) -c $< -o $@;; \
	esac
