# -*- Makefile -*-

VPATH	= $(BINDIR)
exec	= $(BINDIR)/almnfile2rmap
src	= almnfile2rmap.f90
obj 	= almnfile2rmap.o

%.o: %.f90 $(inc)
	$(FC) $(FFLAGS) $(PIXFLAGS) $(LDFLAGS) -c $< -o $@

default: $(exec)

$(exec): $(obj) $(LIBDIR)/lib3dex.a
	$(FC) $(FFLAGS) -o $(exec) $(obj) $(LDFLAGS) $(PIXFLAGS)
