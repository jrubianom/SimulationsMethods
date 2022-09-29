CXX = g++
FLAGS = -std=c++17 -O3
RUN = ./
SOURCES = $(wildcard *.cpp)
DEPENDENCIES = $(SOURCES:%.cpp=.objects/%.o)

SHARE_DIR = /media/sf_ArchData/

.PHONY: main send clean oclean

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

send:
	@echo -e 'Moving graphs ... \c'
	@rm -rf $(SHARE_DIR)/$(notdir $(shell pwd))
	@cp -r data $(SHARE_DIR)/$(notdir $(shell pwd))
	@echo -e 'Done!'

clean:
	@rm -rf data/* *.x

oclean:
	@rm -rf .objects/*.o
