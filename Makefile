SHELL := /bin/bash
CPP := g++

DIRS := lib bin

CFLAGS  := -std=c++11 -Wall -O3 -Isrc
CGFLAGS := -std=c++11 -Wall -g  -Isrc

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --libs)

LHAPDF_CFLAGS := $(shell lhapdf-config --cppflags)
LHAPDF_LIBS   := $(shell lhapdf-config --ldflags)

.PHONY: all clean

EXE := bin/test_rew_calc

all: $(DIRS) $(EXE)

# directories rule
$(DIRS):
	@mkdir -p $@

# class object rules
lib/BHEvent.o lib/SJClusterAlg.o lib/weight_ent.o lib/hist_wrap.o: lib/%.o: src/%.cc src/%.h
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m ... "
	@$(CPP) $(CGFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

lib/rew_calc.o: lib/%.o: src/%.cc src/%.h
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m ... "
	@$(CPP) $(CGFLAGS) $(ROOT_CFLAGS) $(LHAPDF_CFLAGS) -c $(filter %.cc,$^) -o $@

# main object rules
lib/reweigh.o: lib/%.o: src/%.cc
	@echo -e "Compiling \E[0;49;94m"$@"\E[0;0m ... "
	@$(CPP) $(CGFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

lib/test_rew_calc.o: lib/%.o: test/%.cc
	@echo -e "Compiling \E[0;49;94m"$@"\E[0;0m ... "
	@$(CPP) $(CGFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

# executable rules
bin/reweigh: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) $(LHAPDF_LIBS) -lboost_program_options

bin/test_rew_calc: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) $(LHAPDF_LIBS)

# OBJ dependencies
lib/rew_calc.o     : src/BHEvent.h

# EXE_OBJ dependencies
lib/test_rew_calc.o: src/rew_calc.h src/BHEvent.h
#lib/reweigh.o      : src/BHEvent.h src/reweighter.h

# EXE dependencies
bin/test_rew_calc  : lib/rew_calc.o lib/BHEvent.o
#bin/reweigh        : lib/BHEvent.o lib/reweighter.o

clean:
	rm -rf bin lib

