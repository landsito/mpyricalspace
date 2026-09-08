# Regenerate the f2py signature (the .pyf meson consumes) from the IRI-2016 sources
# + iri_batch.f90:
#   make sig      -> iri16f.pyf
#   make mod      -> build the extension in-place (quick iteration; the real build is meson)
#
# The vendor .for files compile as-is -- no patch step.
fortran_files = irisub.for irifun.for cira.for igrf.for iridreg.for iritec.for \
                iriflip.for cosd_sind.for iri_batch.f90
mod  = iri16f
pyf  = iri16f.pyf
only = iri_sub read_ig_rz readapf107 iri_batch

FFLAGS = -std=legacy -w -O2 -fbacktrace -fno-automatic -fPIC
FC     = gfortran

.PHONY: sig mod clean

sig:
	f2py --overwrite-signature -m $(mod) -h $(pyf) $(fortran_files) only: $(only) :

mod: sig
	f2py -c $(pyf) $(fortran_files) --fcompiler=$(FC) --f77flags="$(FFLAGS)"

clean:
	$(RM) *.so *.mod
