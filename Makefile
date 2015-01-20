SHELL := /bin/bash
CPP := g++

DIRS := lib bin

CFLAGS := -std=c++11 -Wall -O3 -Itools/include -Iparts/include

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --libs)

LHAPDF_CFLAGS := $(shell lhapdf-config --cppflags)
LHAPDF_LIBS   := $(shell lhapdf-config --ldflags)

FJ_LIBS := -lfastjet

.PHONY: all misc tools parts clean deepclean

all: $(DIRS) bin/reweigh bin/hist_H2j bin/hist_H3j bin/plot bin/test_H3j bin/test_fj_H3j bin/inspect_bh

misc: $(DIRS) bin/cross_section_hist bin/cross_section_bh bin/hist_weights bin/select_old_weight_hists bin/draw_together

tools parts:
	@$(MAKE) -C $@

# directories rule
lib bin:
	@mkdir -p $@

# main object rules
lib/cross_section_bh.o lib/cross_section_hist.o lib/test_rew_calc.o lib/reweigh.o lib/hist_weights.o lib/select_old_weight_hists.o lib/draw_together.o lib/plot.o lib/test_H3j.o lib/test_fj_H3j.o lib/inspect_bh.o: lib/%.o: src/%.cc
	@echo -e "Compiling \E[0;49;94m"$@"\E[0;0m"
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

lib/hist_H2j.o lib/hist_H3j.o: lib/%.o: src/%.cc
	@echo -e "Compiling \E[0;49;94m"$@"\E[0;0m"
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -DCONFDIR="\"`pwd -P`/config\"" \
		-c $(filter %.cc,$^) -o $@

# executable rules
bin/cross_section_bh bin/cross_section_hist bin/inspect_bh: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS)

bin/reweigh: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) $(LHAPDF_LIBS) -lboost_program_options

bin/hist_weights: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) -Wl,--no-as-needed $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_regex

bin/plot bin/select_old_weight_hists: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_program_options -lboost_regex

bin/draw_together bin/test_H3j: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_program_options

bin/test_fj_H3j: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lm -lfastjettools -lfastjet

bin/hist_H2j: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) -Wl,--no-as-needed $(filter %.o,$^) -o $@ $(ROOT_LIBS) -lboost_program_options -lboost_regex

bin/hist_H3j: bin/%: lib/%.o
	@echo -e "Linking \E[0;49;92m"$@"\E[0;0m"
	@$(CPP) -Wl,--no-as-needed $(filter %.o,$^) -o $@ $(ROOT_LIBS) $(FJ_LIBS) -lboost_program_options -lboost_regex

# EXE_OBJ dependencies
lib/cross_section_bh.o: parts/include/BHEvent.h
lib/inspect_bh.o: parts/include/BHEvent.h
lib/test_rew_calc.o: parts/include/rew_calc.h parts/include/BHEvent.h

lib/reweigh.o: tools/include/timed_counter.h parts/include/rew_calc.h parts/include/BHEvent.h

lib/hist_H2j.o lib/hist_H3j.o lib/test_H3j.o lib/test_fj_H3j.o: tools/include/timed_counter.h parts/include/BHEvent.h parts/include/SJClusterAlg.h
lib/plot.o: tools/include/propmap.h

lib/hist_weights.o lib/hist_H2j.o lib/hist_H3j.o: tools/include/csshists.h

# EXE dependencies
bin/cross_section_bh: parts/lib/BHEvent.o
bin/inspect_bh: parts/lib/BHEvent.o
bin/test_rew_calc: parts/lib/rew_calc.o parts/lib/BHEvent.o

bin/reweigh: tools/lib/timed_counter.o parts/lib/rew_calc.o parts/lib/BHEvent.o

bin/hist_H2j bin/hist_H3j bin/test_H3j bin/test_fj_H3j: tools/lib/timed_counter.o parts/lib/BHEvent.o parts/lib/SJClusterAlg.o

bin/hist_weights bin/hist_H2j bin/hist_H3j: tools/lib/csshists.o

# tools rule
tools/lib/%.o: tools/src/%.cc tools/include/%.h
	@$(MAKE) -C tools lib/$*.o

# parts rule
parts/lib/%.o: parts/src/%.cc parts/include/%.h
	@$(MAKE) -C parts lib/$*.o

clean:
	rm -rf bin/* lib/*

deepclean: clean
	@$(MAKE) -C tools clean
	@$(MAKE) -C parts clean
