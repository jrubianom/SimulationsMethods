CXX = g++
FLAGS = -std=c++11 -O3
RUN = ./
SOURCES = $(wildcard *.cpp)
DEPENDENCIES = $(SOURCES:%.cpp=.objects/%.o)

.PHONY: main clean oclean

main: main.x
	@rm -rf data/*
	@echo -e 'Running program ... \n'
	@$(RUN)$<
	@echo -e '\n\nDone!\n'

main.x: $(DEPENDENCIES)
	@echo -e 'Compiling' $@ '... \c'
	@$(CXX) $(FLAGS) $^ -o $@
	@echo -e 'Done!\n'

.objects/%.o: %.cpp
	@echo -e 'Building' $@ '... \c'
	@$(CXX) $(FLAGS) -c $< -o $@
	@echo -e 'Done!\n'

mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
mkfile_dir := $(dir $(mkfile_path))

mkfile := $(basename $(mkfile_dir))

aux = $(shell pwd)

graph:
	echo $(abspath $(lastword $(MAKEFILE_LIST)))
	echo $(basename $(abspath $(lastword $(MAKEFILE_LIST))))
	echo $(basename /home/3a/)
	echo $(basename $(aux))

clean:
	@rm -rf data/* *.x

oclean:
	@rm -rf .objects/*.o
