# Two products:
#   make            -> eefm_plot         (Patrick Alken's original standalone tool)
#   make pymod      -> eefm1<EXT>.so     (the CPython extension meson installs; handy
#                                         for iterating on eefm1.c without a full build)
CC      = gcc
CFLAGS  = -g -O2 -Wall -I.
LDFLAGS = -lm -lgsl -lgslcblas

MODEL_SRC = eefm_basis.c eefm_calc.c ndlinear.c
OBJECTS   = eefm_basis.o eefm_calc.o eefm_plot.o ndlinear.o

PYINC := $(shell python -c "import sysconfig, numpy; print('-I' + sysconfig.get_path('include'), '-I' + numpy.get_include())")
PYEXT := $(shell python -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")

.PHONY: all pymod clean

all: $(OBJECTS)
	$(CC) -o eefm_plot $(OBJECTS) $(LDFLAGS)

pymod: eefm1.c $(MODEL_SRC)
	$(CC) -shared -fPIC $(CFLAGS) $(PYINC) $^ -o eefm1$(PYEXT) $(LDFLAGS)

.SUFFIXES: .c.o
.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o eefm_plot eefm_plot.exe *.so
