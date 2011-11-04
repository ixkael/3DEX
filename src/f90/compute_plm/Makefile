# -*- Makefile -*-

VPATH	= $(BINDIR)
exec	= $(BINDIR)/compute_plm
src	= compute_plm.f90
obj 	= compute_plm.o

%.o: %.f90 $(inc)
	$(FC) $(FFLAGS) $(PIXFLAGS) $(LDFLAGS) -c $< -o $@

default: $(exec)

$(exec): $(obj) $(LIBDIR)/lib3dex.a
	$(FC) $(FFLAGS) -o $(exec) $(obj) $(LDFLAGS) $(PIXFLAGS)
