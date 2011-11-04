# -*- Makefile -*-

VPATH	= $(BINDIR)
exec	= $(BINDIR)/compute_qln
src	= compute_qln.f90
obj 	= compute_qln.o

%.o: %.f90 $(inc)
	$(FC) $(FFLAGS) $(PIXFLAGS) $(LDFLAGS) -c $< -o $@

default: $(exec)

$(exec): $(obj) $(LIBDIR)/lib3dex.a
	$(FC) $(FFLAGS) -o $(exec) $(obj) $(LDFLAGS) $(PIXFLAGS)
