# -*- Makefile -*-

obj0 = ../external/quadpack.o ../external/quadrule.o
obj1 = ../external/rfsfns.o
obj2 = f3dex_stats.o f3dex_utils.o f3dex_fitstools.o
obj3 = f3dex_cosmotools.o f3dex_transforms.o

obj   = $(obj0) $(obj1) $(obj2) $(obj3)

default: lib3dex.a 

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

f3dex_cosmotools.o: $(obj2)
f3dex_utils.o: $(obj1)
f3dex_stats.o: $(obj1)
f3dex_transforms.o: $(obj2)
f3dex_fitstools.o: $(obj1)
../external/rfsfns.o: $(obj0)

lib3dex.a: $(obj)
	$(AR) lib3dex.a $(obj)
	@cp -p lib3dex.a $(LIBDIR)/
	@for p in *.$(MOD); do \
	cp $$p $(INCDIR)/; \
	done

