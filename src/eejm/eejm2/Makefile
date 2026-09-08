# Two products:
#   make            -> eej_plot          (Patrick Alken's original standalone tool)
#   make pymod      -> eejm2<EXT>.so     (the CPython extension meson installs; handy
#                                         for iterating on eejm2.c without a full build)
CC      = gcc
CFLAGS  = -g -O2 -Wall -I.
LDFLAGS = -lm -lgsl -lgslcblas

MODEL_SRC = eej_basis.c eej_calc.c ndlinear.c common.c coord.c
OBJECTS   = eej_basis.o eej_calc.o eej_plot.o ndlinear.o common.o coord.o

PYINC := $(shell python -c "import sysconfig, numpy; print('-I' + sysconfig.get_path('include'), '-I' + numpy.get_include())")
PYEXT := $(shell python -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")

.PHONY: all pymod clean

all: $(OBJECTS)
	$(CC) -o eej_plot $(OBJECTS) $(LDFLAGS)

pymod: eejm2.c $(MODEL_SRC)
	$(CC) -shared -fPIC $(CFLAGS) $(PYINC) $^ -o eejm2$(PYEXT) $(LDFLAGS)

.SUFFIXES: .c.o
.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o eej_plot eej_plot.exe *.so
