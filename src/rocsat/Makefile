fortran_files = vfjmodroc1.F90 
mod = rocsatf
only = getverticaldrift 

FFLAGS = -std=legacy -ffixed-line-length-none -undefined dynamic_lookup
FC=gfortran

.PHONY: clean all

clean:
	$(RM) *.pyf *.so

%.pyf:
	f2py --overwrite-signature -m $(mod) -h $@ $(fortran_files) only: $(only) :

$(mod): $(mod).pyf
	f2py -c $< $(fortran_files) --fcompiler=$(FC) --f77flags="$(FFLAGS)"

all: clean $(mod)