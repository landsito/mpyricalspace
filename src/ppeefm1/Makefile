# Two products:
#   make            -> rteef             (Patrick Alken's original standalone tool)
#   make pymod      -> ppeefm1<EXT>.so   (the CPython extension meson installs; handy
#                                         for iterating on ppeefm1.c without a full build)
CC      = gcc
CFLAGS  = -g -O2 -Wall -I. -Indlinear
LDFLAGS = -lm -lgsl -lgslcblas

MODEL_SRC = rteef.c eefm_basis.c eefm_calc.c ndlinear.c f107.c common.c lunar.c
OBJECTS   = rteef.o eefm_basis.o eefm_calc.o ndlinear.o f107.o common.o lunar.o

PYINC := $(shell python -c "import sysconfig, numpy; print('-I' + sysconfig.get_path('include'), '-I' + numpy.get_include())")
PYEXT := $(shell python -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")

.PHONY: all pymod clean

all: $(OBJECTS)
	$(CC) -o rteef $(OBJECTS) $(LDFLAGS)

pymod: ppeefm1.c $(MODEL_SRC)
	$(CC) -shared -fPIC $(CFLAGS) $(PYINC) $^ -o ppeefm1$(PYEXT) $(LDFLAGS)

.SUFFIXES: .c.o
.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o rteef *.so
