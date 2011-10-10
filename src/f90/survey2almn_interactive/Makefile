# -*- Makefile -*-

VPATH	= $(BINDIR)
exec	= $(BINDIR)/survey2almn_interactive
src	= survey2almn_interactive.f90
obj 	= survey2almn_interactive.o

%.o: %.f90 $(inc)
	$(FC) $(FFLAGS) $(PIXFLAGS) $(LDFLAGS) -c $< -o $@

default: $(exec)

$(exec): $(obj) $(LIBDIR)/lib3dex.a 
	$(FC) $(FFLAGS) -o $(exec) $(obj) $(LDFLAGS) $(PIXFLAGS)
