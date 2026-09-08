# Regenerate the f2py signature (the .pyf meson consumes) from igrf11.f:
#   make sig      -> igrf11.pyf
#   make mod      -> build the extension in-place (quick iteration; the real build is meson)
fortran_files = igrf11.f
mod  = igrf11f
pyf  = igrf11.pyf
only = igrf11syn

FFLAGS = -std=legacy
FC     = gfortran

.PHONY: sig mod clean

sig:
	f2py --overwrite-signature -m $(mod) -h $(pyf) $(fortran_files) only: $(only) :

mod: sig
	f2py -c $(pyf) $(fortran_files) --fcompiler=$(FC) --f77flags="$(FFLAGS)"

clean:
	$(RM) *.so *.mod
