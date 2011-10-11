# -*- Makefile -*-

VPATH	= $(BINDIR)
exec	= $(BINDIR)/almn2cln
src	= almn2cln.f90
obj 	= almn2cln.o

%.o: %.f90 $(inc)
	$(FC) $(FFLAGS) $(PIXFLAGS) $(LDFLAGS) -c $< -o $@

default: $(exec)

$(exec): $(obj) $(LIBDIR)/lib3dex.a
	$(FC) $(FFLAGS) -o $(exec) $(obj) $(LDFLAGS) $(PIXFLAGS)
