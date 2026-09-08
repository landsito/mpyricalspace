# Regenerate the f2py signature (the .pyf meson consumes) from the HWM07 sources
# + hwm07_batch.f90:
#   make sig      -> hwm07f.pyf
#   make mod      -> build the extension in-place (quick iteration; the real build is meson)
#
# dwm07b.f90 / apexcord.f90 are the NRL release, unchanged.  hwm07e.f90 has one fix:
# subroutine HWMupdate indexes ebz/ebm directly instead of via pointer-aliased
# `target` dummies (non-conforming F95).  See the file header.
fortran_files = hwm07e.f90 dwm07b.f90 apexcord.f90 hwm07_batch.f90
mod  = hwm07f
pyf  = hwm07f.pyf
only = hwm07 dwm07b_hwm_interface hwm07_batch

FFLAGS = -std=legacy -fno-automatic
FC     = gfortran

.PHONY: sig mod clean

sig:
	f2py --overwrite-signature -m $(mod) -h $(pyf) $(fortran_files) only: $(only) :

mod: sig
	f2py -c $(pyf) $(fortran_files) --fcompiler=$(FC) --f90flags="$(FFLAGS)"

clean:
	$(RM) *.so *.mod
