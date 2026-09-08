# Regenerate the f2py signature (the .pyf meson consumes) from the IRI-2012 sources
# + iri_batch.f90:
#   make sig      -> iri12f.pyf
#   make mod      -> build the extension in-place (quick iteration; the real build is meson)
#
# The vendor .for files compile as-is (no patches). IRI-2012 reads its index files
# (ig_rz.dat / apf107.dat) on the first call, so read_ig_rz / readapf107 are not exported.
fortran_files = irisub.for cira.for igrf.for iridreg.for irifun.for iritec.for \
                iriflip.for iri_batch.f90
mod  = iri12f
pyf  = iri12f.pyf
only = iri_sub iri_batch

FFLAGS = -std=legacy -w -O2 -fbacktrace -fno-automatic -fPIC
FC     = gfortran

.PHONY: sig mod clean

sig:
	f2py --overwrite-signature -m $(mod) -h $(pyf) $(fortran_files) only: $(only) :

mod: sig
	f2py -c $(pyf) $(fortran_files) --fcompiler=$(FC) --f77flags="$(FFLAGS)"

clean:
	$(RM) *.so *.mod
