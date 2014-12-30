SHELL := /bin/bash
CPP := g++

DIRS := lib bin

CFLAGS  := -std=c++11 -Wall -O3 -Isrc

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --libs)

LHAPDF_CFLAGS := $(shell lhapdf-config --cppflags)
LHAPDF_LIBS   := $(shell lhapdf-config --ldflags)

.PHONY: all misc clean

all: $(DIRS) bin/reweigh bin/hist_H2j bin/plot

misc: $(DIRS) bin/cross_section bin/hist_weights bin/select_old_weight_hists bin/draw_together

# directories rule
$(DIRS):
	@mkdir -p $@

# class object rules
lib/BHEvent.o lib/SJClusterAlg.o lib/weight.o lib/selector_base.o: lib/%.o: src/%.cc src/%.h
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m ... "
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

lib/rew_calc.o: lib/%.o: src/%.cc src/%.h
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m ... "
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) $(LHAPDF_CFLAGS) -c $(filter %.cc,$^) -o $@

# analysis plugins rules
lib/selector_H2j.o: lib/%.o: src/%.cc src/selector.h
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m ... "
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

# main object rules
lib/cross_section_bh.o lib/cross_section_hist.o lib/test_rew_calc.o lib/reweigh.o lib/hist_weights.o lib/select_old_weight_hists.o lib/draw_together.o lib/plot.o: lib/%.o: src/%.cc
	@echo -e "Compiling \E[0;49;94m"$@"\E[0;0m ... "
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

lib/hist.o: lib/%.o: src/%.cc
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

bin/draw_together: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m ... "
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_program_options

bin/hist_H2j: bin/hist_%: lib/hist.o lib/selector_%.o
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
lib/hist.o         : src/BHEvent.h src/SJClusterAlg.h src/weight.h src/timed_counter.h
lib/selector_base.o: src/BHEvent.h src/SJClusterAlg.h src/weight.h src/timed_counter.h
lib/plot.o         : src/propmap11.h

# EXE dependencies
bin/cross_section_bh: lib/BHEvent.o
bin/test_rew_calc  : lib/rew_calc.o lib/BHEvent.o
bin/reweigh        : lib/rew_calc.o lib/BHEvent.o
bin/hist_H2j       : lib/BHEvent.o lib/SJClusterAlg.o lib/weight.o lib/selector_base.o

clean:
	rm -rf bin lib
