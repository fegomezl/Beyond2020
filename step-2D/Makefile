include ../local_config.mk

#Mesh parmeters
SCRIPT=$(shell sed -n 2p data/parameters.txt | cut -d '#' -f 1)
RMIN=$(shell sed -n 3p data/parameters.txt | tr -d -c 0-9.)
RMAX=$(shell sed -n 4p data/parameters.txt | tr -d -c 0-9.)
ZMIN=$(shell sed -n 5p data/parameters.txt | tr -d -c 0-9.)
ZMAX=$(shell sed -n 6p data/parameters.txt | tr -d -c 0-9.)

REF = $(shell sed -n 11p data/parameters.txt | tr -d -c 0-9.)
ORDER = $(shell sed -n 12p data/parameters.txt | tr -d -c 0-9.)
T_FU = $(shell sed -n 15p data/parameters.txt | tr -d -c 0-9.)
ALPHA_L = $(shell sed -n 16p data/parameters.txt | tr -d -c 0-9.)
ALPHA_S = $(shell sed -n 17p data/parameters.txt | tr -d -c 0-9.)
DT = $(shell sed -n 18p data/parameters.txt | tr -d -c 0-9.)
T_FI = $(shell sed -n 19p data/parameters.txt | tr -d -c 0-9.)
VIS = $(shell sed -n 20p data/parameters.txt | tr -d -c 0-9.)
ODE = $(shell sed -n 21p data/parameters.txt | tr -d -c 0-9.)
RELT = $(shell sed -n 22p data/parameters.txt | tr -d -c 0-9.)
ABST = $(shell sed -n 23p data/parameters.txt | tr -d -c 0-9.)


CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS)
RUN = mpirun -np $(PROCCESORS) ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)
MESH = graph/mesh.msh

.PHONY: all graph clean oclean

all: main.x
	@echo -e 'Reading parameters ... \c'
	@bash scripts/parameters_script.sh
	@echo -e 'Generating mesh ... \n'
	@gmsh $(SCRIPT) -format msh2 -o $(MESH) -3 > /dev/null
	@echo -e '\nDone!\n'
	@echo -e 'Running program ... \n'
	@$(RUN)$< --mesh $(MESH) -Rmin $(RMIN) -Rmax $(RMAX) -Zmin $(ZMIN)  -Zmax $(ZMAX) \
						-o $(ORDER) -r $(REF) -T_f $(T_FU) -a_l $(ALPHA_L) -a_s $(ALPHA_S) \
						-dt $(DT) -t_f $(T_FI) -v_s $(VIS) -ode $(ODE) -reltol $(RELT) -abstol $(ABST)
	@echo -e '\nDone!\n'

graph:
	@echo -e 'Moving graphs ... \c'
	@rm -rf $(SHARE_DIR)/step-2D
	@cp -r graph/step-2D $(SHARE_DIR)/step-2D
	@echo -e 'Done!'

main.x: $(DEPENDENCIES)
	@echo -e 'Compiling' $@ '... \c'
	@$(CXX) $(FLAGS) $^ $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

.objects/%.o: code/%.cpp
	@echo -e 'Building' $@ '... \c'
	@$(CXX) $(FLAGS) -c $< $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

clean:
	@rm -rf *.x graph/*

oclean:
	@rm -rf .objects/*.o
