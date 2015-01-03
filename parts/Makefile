SHELL := /bin/bash
CPP := g++

CFLAGS := -std=c++11 -Wall -O3 -Iinclude

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --libs)

LHAPDF_CFLAGS := $(shell lhapdf-config --cppflags)
LHAPDF_LIBS   := $(shell lhapdf-config --ldflags)

.PHONY: all clean

all: lib bin \
	   lib/BHEvent.o lib/rew_calc.o lib/SJClusterAlg.o \
		 bin/test_rew_calc

# directories rule
lib bin:
	@mkdir -p $@

# class object rules
lib/BHEvent.o lib/SJClusterAlg.o: lib/%.o: src/%.cc include/%.h
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m"
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

lib/rew_calc.o: lib/%.o: src/%.cc include/%.h
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m"
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) $(LHAPDF_CFLAGS) -c $(filter %.cc,$^) -o $@

# text executable rules
bin/test_rew_calc: bin/%: test/%.cc
	@echo -e "Compiling \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) $(filter %.cc,$^) -o $@ \
	 -Wl,--no-as-needed $(ROOT_LIBS) $(LHAPDF_LIBS) $(filter %.o,$^)

# dependencies
lib/rew_calc.o     : include/BHEvent.h

bin/test_rew_calc  : include/BHEvent.h lib/rew_calc.o lib/BHEvent.o

clean:
	rm -rf lib/* bin/*
