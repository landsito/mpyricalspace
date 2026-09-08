# Regenerate the f2py signature (the .pyf meson consumes) from hwm93.f + hwm93_batch.f90:
#   make sig      -> hwm93f.pyf
#   make mod      -> build the extension in-place (quick iteration; the real build is meson)
#
# hwm93.f is the NRL release with its non-Fortran header/footer text removed and the
# ISDATE/ISTIME/NAME DATA statements switched to 4H Hollerith -- no model-logic changes.
fortran_files = hwm93.f hwm93_batch.f90
mod  = hwm93f
pyf  = hwm93f.pyf
only = gws5 gws5_batch

FFLAGS = -std=legacy -fno-automatic
FC     = gfortran

.PHONY: sig mod clean

sig:
	f2py --overwrite-signature -m $(mod) -h $(pyf) $(fortran_files) only: $(only) :

mod: sig
	f2py -c $(pyf) $(fortran_files) --fcompiler=$(FC) --f77flags="$(FFLAGS)"

clean:
	$(RM) *.so *.mod
