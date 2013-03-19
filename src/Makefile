CCC = g++

ifeq ($(CBLAS),)
CCOPT = -I/usr/share/R/include -fPIC
GSL = -lrt
else
CCOPT = -I/usr/share/R/include -fPIC -DUSE_CBLAS
GSL = -lgslcblas -lrt
endif

CCFLAGS = $(CCOPT)

ifeq ($(DEB),0)
  DFLAGS = -DNDEBUG -O3 -DOPTIMIZE
else
ifeq ($(DEB),)
  DFLAGS = -DNDEBUG -O3 -DOPTIMIZE
else
ifeq ($(DEB),1)
  DFLAGS = -O3 -Wall -ansi -DCPURE
else
  DFLAGS = -ggdb -g -g3 -DDEB=$(DEB) -Wall -ansi -DCPURE
endif
endif
endif

ifeq ($(CPURE),)
  BIN = ../bin/nmf.so
  DDFLAGS = $(DFLAGS) -L/usr/lib64/R/lib -lR -shared
else
  BIN = ../bin/nmf-standalone
  DDFLAGS = $(DFLAGS) -DCPURE
endif

cpps :=  $(wildcard *.cc)
objects := $(foreach file, $(cpps), $(basename $(file)).o)

all: $(BIN)

obj: $(objects)

$(BIN): $(objects)
	$(CCC) $(CCFLAGS) $(DDFLAGS) -o $(BIN) $(objects) $(GSL)

$(objects): %.o: %.cc %.hh
	$(CCC) $(CCFLAGS) $(DDFLAGS) -c $< -o $@

.PHONY: clean cleanbin realclean cleandist distclean
cleandist: realclean
distclean: realclean
realclean: cleanbin clean

cleanbin:
	-rm -f ../bin/nmf.so ../bin/nmf-standalone
clean:
	-rm -f *.*~ $(objects)
