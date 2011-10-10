# -*- Makefile -*-

obj0 = quadpack.o 
obj1 = rfsfns.o
obj2 = cosmotools.o f3dex.o

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

cosmotools.o: $(obj0)

f3dex.o: $(obj1)

lib3dex.a: $(obj)
	$(AR) lib3dex.a $(obj)
	@cp -p lib3dex.a $(LIBDIR)/
	@for p in *.$(MOD); do \
	cp $$p $(INCDIR)/; \
	done

