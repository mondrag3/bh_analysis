SHELL := /bin/bash
CPP := g++

DIRS := lib bin

CFLAGS  := -std=c++11 -Wall -O3 -Isrc

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --libs)

LHAPDF_CFLAGS := $(shell lhapdf-config --cppflags)
LHAPDF_LIBS   := $(shell lhapdf-config --ldflags)

.PHONY: all misc clean

# bin/test_rew_calc
EXE := bin/reweigh bin/hist_H2j bin/plot bin/test_H3j bin/test_fj_H3j

all: $(DIRS) $(EXE)

misc: $(DIRS) bin/cross_section bin/hist_weights bin/select_old_weight_hists bin/draw_together

# directories rule
$(DIRS):
	@mkdir -p $@

# class object rules
lib/BHEvent.o lib/SJClusterAlg.o: lib/%.o: src/%.cc src/%.h
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m ... "
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

lib/rew_calc.o: lib/%.o: src/%.cc src/%.h
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m ... "
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) $(LHAPDF_CFLAGS) -c $(filter %.cc,$^) -o $@

# main object rules
lib/cross_section_bh.o lib/cross_section_hist.o lib/test_rew_calc.o lib/reweigh.o lib/hist_weights.o lib/select_old_weight_hists.o lib/draw_together.o lib/plot.o lib/test_H3j.o lib/test_fj_H3j.o: lib/%.o: src/%.cc
	@echo -e "Compiling \E[0;49;94m"$@"\E[0;0m ... "
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

lib/hist_H2j.o lib/hist_H3j.o: lib/%.o: src/%.cc
	@echo -e "Compiling \E[0;49;94m"$@"\E[0;0m ... "
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -DCONFDIR="\"`pwd -P`/config\"" \
		-c $(filter %.cc,$^) -o $@

# executable rules
bin/cross_section_bh bin/cross_section_hist: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS)

bin/reweigh: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) $(LHAPDF_LIBS) -lboost_program_options

bin/test_rew_calc: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) $(LHAPDF_LIBS)

bin/hist_weights: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) -Wl,--no-as-needed $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_regex -lkiwihist

bin/select_old_weight_hists: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_program_options -lboost_regex

bin/draw_together bin/test_H3j: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_program_options

bin/test_fj_H3j: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lm -lfastjettools -lfastjet

bin/hist_H2j bin/hist_H3j: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) -Wl,--no-as-needed $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_program_options -lboost_regex -lkiwihist

bin/plot: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_program_options -lboost_regex

# OBJ dependencies
lib/rew_calc.o     : src/BHEvent.h

# EXE_OBJ dependencies
lib/cross_section_bh.o: src/BHEvent.h
lib/test_rew_calc.o: src/rew_calc.h src/BHEvent.h
lib/reweigh.o      : src/rew_calc.h src/BHEvent.h src/timed_counter.h
lib/hist_H2j.o lib/hist_H3j.o : src/BHEvent.h src/SJClusterAlg.h src/finder.h src/timed_counter.h
lib/test_H3j.o lib/test_fj_H3j.o : src/BHEvent.h src/SJClusterAlg.h src/timed_counter.h
lib/plot.o         : src/propmap11.h

# EXE dependencies
bin/cross_section_bh: lib/BHEvent.o
bin/test_rew_calc  : lib/rew_calc.o lib/BHEvent.o
bin/reweigh        : lib/rew_calc.o lib/BHEvent.o
bin/hist_H2j bin/hist_H3j bin/test_H3j bin/test_fj_H3j : lib/BHEvent.o lib/SJClusterAlg.o

clean:
	rm -rf bin lib
