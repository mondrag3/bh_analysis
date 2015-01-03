SHELL := /bin/bash
CPP := g++

CFLAGS := -std=c++11 -Wall -O3 -Iinclude

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --libs)

.PHONY: all clean

all: lib lib/timed_counter.o lib/csshists.o

# directories rule
lib bin:
	@mkdir -p $@

# class object rules
lib/timed_counter.o: lib/%.o: src/%.cc include/%.h
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m"
	@$(CPP) $(CFLAGS) -c $(filter %.cc,$^) -o $@

lib/csshists.o: lib/%.o: src/%.cc include/%.h
	@echo -e "Compiling \E[0;49;96m"$@"\E[0;0m"
	@$(CPP) $(CFLAGS) $(ROOT_CFLAGS) -c $(filter %.cc,$^) -o $@

clean:
	rm -rf lib/* bin/*
