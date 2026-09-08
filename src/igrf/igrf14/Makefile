# Regenerate the f2py signature (the .pyf meson consumes) from igrf14.f:
#   make sig      -> igrf14.pyf
#   make mod      -> build the extension in-place (quick iteration; the real build is meson)
fortran_files = igrf14.f
mod  = igrf14f
pyf  = igrf14.pyf
only = igrf14syn

FFLAGS = -std=legacy
FC     = gfortran

.PHONY: sig mod clean

sig:
	f2py --overwrite-signature -m $(mod) -h $(pyf) $(fortran_files) only: $(only) :

mod: sig
	f2py -c $(pyf) $(fortran_files) --fcompiler=$(FC) --f77flags="$(FFLAGS)"

clean:
	$(RM) *.so *.mod
