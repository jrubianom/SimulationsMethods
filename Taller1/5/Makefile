#Compiling parameters
PDF_VIEWER=okular
FLAGS = -std=c++17
SOURCES = $(shell find . -type f -name "*.cpp" | sed -e 's/\.cpp//g' -e 's/\.\///g')
BCyan=\033[1;36m
NC=\033[0m

.PHONY: graph
.PRECIOUS: %.x

all:
	@echo -e "${BCyan}PROGRAMAS:${NC}"
	@echo $(SOURCES) | tr " " "\n"
	@echo -e "${BCyan}GRAPH:${NC}"
	@echo "make graph gp=gnuplot_script corre el script gp.gp y abre el PDF resultado gp.pdf con PDF_VIEWER (okular)"
	@echo -e "${BCyan}CLEAN:${NC}"
	@echo "make clean borra los .x y .dat en todas las subcarpetas"

punto_5a: punto_5a.cpp punto_5a.gp
	@echo 'Compiling $@.cpp'
	@g++ $(FLAGS) $< -o $@.x
	@echo -e 'Running $< \nData stored in data/$@.dat'
	@./$@.x > data/$@.dat

punto_5b: punto_5b.cpp punto_5b.gp
	@echo 'Compiling $@.cpp'
	@g++ $(FLAGS) $< -o $@.x
	@echo -e 'Running $< \nData stored in data/$@.dat'
	@./$@.x > data/$@.dat

punto_5f: punto_5f.cpp punto_5f.gp
	@echo 'Compiling $@.cpp'
	@g++ $(FLAGS)  -lgsl -lgslcblas -lm $< -o $@.x
	@echo -e 'Running $< \nData stored in data/$@.dat'
	@./$@.x > data/$@.dat

punto_5f_histogram: punto_5f_histogram.cpp punto_5f_histogram.gp
	@echo 'Compiling $@.cpp'
	@g++ $(FLAGS) -lgsl -lgslcblas -lm $< -o $@.x
	@echo -e 'Running $< \nData stored in data/$@.dat'
	@./$@.x > data/$@.dat

punto_5h: punto_5h.cpp punto_5h.gp
	@echo 'Compiling $@.cpp'
	@g++ $(FLAGS) -lgsl -lgslcblas -lm $< -o $@.x
	@echo -e 'Running $< \nData stored in data/$@.dat'
	@./$@.x > data/$@.dat

punto_5i: punto_5i.cpp punto_5i.gp
	@echo 'Compiling $@.cpp'
	@g++ $(FLAGS) -lgsl -lgslcblas -lm $< -o $@.x
	@echo -e 'Running $< \nData stored in data/$@.dat'
	@./$@.x > data/$@.dat



clean:
	@find . -type f -name "*.x" -delete
	@find . -type f -name "*.dat" -delete
	@find . -type f -name "*.log" -delete
