# -*- Makefile -*-

VPATH	= $(BINDIR)
exec	= $(BINDIR)/almnfile2cln
src	= almnfile2cln.f90
obj 	= almnfile2cln.o

%.o: %.f90 $(inc)
	$(FC) $(FFLAGS) $(PIXFLAGS) $(LDFLAGS) -c $< -o $@

default: $(exec)

$(exec): $(obj) $(LIBDIR)/lib3dex.a
	$(FC) $(FFLAGS) -o $(exec) $(obj) $(LDFLAGS) $(PIXFLAGS)
