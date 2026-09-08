# Regenerate scherliessfejerf.pyf from the !f2py directives in StormDrift.F90:
#   make pyf
fortran_files = StormDrift.F90
mod  = scherliessfejerf
only = vdrift_quiet_model vdrift_storm_model

FFLAGS = -std=legacy
FC     = gfortran

.PHONY: pyf mod clean

pyf:
	f2py --overwrite-signature -m $(mod) -h $(mod).pyf $(fortran_files) only: $(only) :

mod: pyf
	f2py -c $(mod).pyf $(fortran_files) --fcompiler=$(FC) --f77flags="$(FFLAGS)"

clean:
	$(RM) $(mod).pyf *.so *.mod
