SHELL := /bin/bash
CPP := g++

DIRS := lib bin

CFLAGS  := -std=c++11 -Wall -O3 -Isrc

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --libs)

LHAPDF_CFLAGS := $(shell lhapdf-config --cppflags)
LHAPDF_LIBS   := $(shell lhapdf-config --ldflags)

.PHONY: all clean

# bin/test_rew_calc
EXE := bin/reweigh \
       bin/hist_weights bin/select_old_weight_hists bin/draw_together \
       bin/read_withFriending

all: $(DIRS) $(EXE)

# directories rule
$(DIRS):
	@mkdir -p $@

# class object rules
lib/BHEvent.o lib/SJClusterAlg.o lib/weight_ent.o lib/hist_wrap.o: lib/%.o: src/%.cc src/%.h
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m ... "
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

lib/rew_calc.o: lib/%.o: src/%.cc src/%.h
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m ... "
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) $(LHAPDF_CFLAGS) -c $(filter %.cc,$^) -o $@

# main object rules
lib/test_rew_calc.o lib/reweigh.o lib/hist_weights.o lib/select_old_weight_hists.o lib/draw_together.o lib/read_withFriending.o: lib/%.o: src/%.cc
	@echo -e "Compiling \E[0;49;94m"$@"\E[0;0m ... "
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

# executable rules
bin/reweigh: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) $(LHAPDF_LIBS) -lboost_program_options

bin/test_rew_calc: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) $(LHAPDF_LIBS)

bin/hist_weights: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_regex

bin/select_old_weight_hists: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_program_options -lboost_regex

bin/draw_together: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_program_options

bin/read_withFriending: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) -Wl,--no-as-needed $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_program_options -lboost_regex -lkiwihist

# OBJ dependencies
lib/rew_calc.o     : src/BHEvent.h

# EXE_OBJ dependencies
lib/test_rew_calc.o: src/rew_calc.h src/BHEvent.h
lib/reweigh.o      : src/rew_calc.h src/BHEvent.h
lib/hist_weights.o : src/hist_wrap.h
lib/select_old_weight_hists.o: src/propmap.h
lib/read_withFriending.o: src/BHEvent.h src/SJClusterAlg.h

# EXE dependencies
bin/test_rew_calc  : lib/rew_calc.o lib/BHEvent.o
bin/reweigh        : lib/rew_calc.o lib/BHEvent.o
bin/hist_weights   : lib/hist_wrap.o
bin/read_withFriending: lib/BHEvent.o lib/SJClusterAlg.o

clean:
	rm -rf bin lib
