# Regenerate the f2py signature (the .pyf meson consumes) from igrf13.f:
#   make sig      -> igrf13.pyf
#   make mod      -> build the extension in-place (quick iteration; the real build is meson)
fortran_files = igrf13.f
mod  = igrf13f
pyf  = igrf13.pyf
only = igrf13syn

FFLAGS = -std=legacy
FC     = gfortran

.PHONY: sig mod clean

sig:
	f2py --overwrite-signature -m $(mod) -h $(pyf) $(fortran_files) only: $(only) :

mod: sig
	f2py -c $(pyf) $(fortran_files) --fcompiler=$(FC) --f77flags="$(FFLAGS)"

clean:
	$(RM) *.so *.mod
