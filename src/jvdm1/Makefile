# Two products:
#   make            -> drift_plot        (Patrick Alken's original standalone tool)
#   make pymod      -> jvdm1<EXT>.so     (the CPython extension meson installs; handy
#                                         for iterating on jvdm1.c without a full build)
CC      = gcc
CFLAGS  = -g -O2 -Wall -I. -Indlinear
LDFLAGS = -lm -lgsl -lgslcblas

MODEL_SRC = drift_basis.c drift_calc.c ndlinear.c
OBJECTS   = drift_basis.o drift_calc.o drift_plot.o ndlinear.o

PYINC := $(shell python -c "import sysconfig, numpy; print('-I' + sysconfig.get_path('include'), '-I' + numpy.get_include())")
PYEXT := $(shell python -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")

.PHONY: all pymod clean

all: $(OBJECTS)
	$(CC) -o drift_plot $(OBJECTS) $(LDFLAGS)

pymod: jvdm1.c $(MODEL_SRC)
	$(CC) -shared -fPIC $(CFLAGS) $(PYINC) $^ -o jvdm1$(PYEXT) $(LDFLAGS)

.SUFFIXES: .c.o
.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o drift_plot drift_plot.exe *.so
