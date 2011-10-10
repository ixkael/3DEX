# -*- Makefile -*-

VPATH	= $(BINDIR)
exec	= $(BINDIR)/survey2almn
src	= survey2almn.f90
obj 	= survey2almn.o

%.o: %.f90 $(inc)
	$(FC) $(FFLAGS) $(PIXFLAGS) $(LDFLAGS) -c $< -o $@

default: $(exec)

$(exec): $(obj) $(LIBDIR)/lib3dex.a
	$(FC) $(FFLAGS) -o $(exec) $(obj) $(LDFLAGS) $(PIXFLAGS)
