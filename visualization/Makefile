include ../local_config.mk

#Mesh parameters
NMESH = $(shell sed -n 2p data/parameters.txt)  
MESH = $(filter-out #Mesh_file,$(NMESH))
OUTER_RADIUS = $(shell sed -n 4p data/parameters.txt | tr -d -c 0-9.)
INNER_RADIUS = $(shell sed -n 5p data/parameters.txt | tr -d -c 0-9.)

#Program parameters
REF = $(shell sed -n 11p data/parameters.txt | tr -d -c 0-9.) 
ORDER = $(shell sed -n 12p data/parameters.txt | tr -d -c 0-9.)
T_FU = $(shell sed -n 15p data/parameters.txt | tr -d -c 0-9.)
T_I = $(shell sed -n 16p data/parameters.txt | tr -d -c 0-9.)
ALPHA_L = $(shell sed -n 17p data/parameters.txt | tr -d -c 0-9.)
ALPHA_S = $(shell sed -n 18p data/parameters.txt | tr -d -c 0-9.)
LAMBDA = $(shell sed -n 19p data/parameters.txt | tr -d -c 0-9.)
DT = $(shell sed -n 20p data/parameters.txt | tr -d -c 0-9.)
T_IN = $(shell sed -n 21p data/parameters.txt | tr -d -c 0-9.)
T_FI = $(shell sed -n 22p data/parameters.txt | tr -d -c 0-9.)

CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS)
RUN = mpirun -np $(PROCCESORS) ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all graph clean oclean

all: main.x data/parameters.txt
	@echo -e 'Reading parameters ... \c'
	@bash script/parameters_script.sh
	@echo -e 'Generating mesh ... \n'
	@gmsh script/mesh_script.geo -format msh2 -o $(MESH) -3
	@echo -e '\nDone!\n'
	@echo -e 'Running program ... \n'
	@$(RUN)$< --mesh $(MESH) \
		 				--temperature_fusion $(T_FU) --temperature_initial $(T_I) \
						--alpha_liquid $(ALPHA_L) --alpha_solid $(ALPHA_S) --lambda $(LAMBDA) \
						--order $(ORDER) --refinements $(REF) --time_step $(DT) \
						--t_init $(T_IN) --t_final $(T_FI)
	@echo -e '\nDone!\n'

graph:
	@echo -e 'Moving graphs ... \c'
	@cp -r graph/ $(SHARE_DIR)
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
