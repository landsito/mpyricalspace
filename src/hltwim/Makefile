# Regenerate hltwimf.pyf from the !f2py directives in hltwim.f90.
#   make pyf
# f2py emits the `use`d modules (base, base1, alf, gd2qdc) as `module` blocks;
# their _f2pyinit<mod>_ symbols don't resolve in meson's flat-namespace build
# (same issue as iri12), so we strip every inner module block and keep only the
# two subroutine wrappers.

fortran_files = hltwim.f90
mod   = hltwimf
only  = hltwim hltwim_batch

.PHONY: pyf clean

pyf:
	f2py --overwrite-signature -m $(mod) -h $(mod).pyf.raw $(fortran_files) only: $(only) :
	awk 'BEGIN{s=0} /^[ \t]+module [a-z]/{s=1} s==0{print} /^[ \t]+end module /{s=0}' \
	    $(mod).pyf.raw > $(mod).pyf
	$(RM) $(mod).pyf.raw

clean:
	$(RM) $(mod).pyf.raw *.so *.mod
