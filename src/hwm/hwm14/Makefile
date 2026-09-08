# Regenerate the f2py signature (the .pyf meson consumes) from hwm14.f90 + hwm14_batch.f90:
#   make sig      -> hwm14f.pyf
#   make mod      -> build the extension in-place (quick iteration; the real build is meson)
fortran_files = hwm14.f90 hwm14_batch.f90
mod  = hwm14f
pyf  = hwm14f.pyf
only = hwm14 dwm07 hwm14_batch

FFLAGS = -std=legacy
FC     = gfortran

.PHONY: sig mod clean

sig:
	f2py --overwrite-signature -m $(mod) -h $(pyf) $(fortran_files) only: $(only) :

mod: sig
	f2py -c $(pyf) $(fortran_files) --fcompiler=$(FC) --f77flags="$(FFLAGS)"

clean:
	$(RM) *.so *.mod
