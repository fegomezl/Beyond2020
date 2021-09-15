include ../../local_config.mk

#Mesh parameters
SCRIPT=$(shell sed -n 2p settings/parameters.txt | cut -d '#' -f 1)
INNER_RADIUS=$(shell sed -n 3p settings/parameters.txt | tr -d -c 0-9.)
OUTER_RADIUS=$(shell sed -n 4p settings/parameters.txt | tr -d -c 0-9.)
HEIGHT=$(shell sed -n 6p settings/parameters.txt | tr -d -c 0-9.)  

#Program parameters
REF=$(shell sed -n 11p settings/parameters.txt | tr -d -c 0-9.) 
ORDER=$(shell sed -n 12p settings/parameters.txt | tr -d -c 0-9.)

#Physical parameters
T_FU=$(shell sed -n 15p settings/parameters.txt | tr -d -c 0-9.-)
DELTA_T=$(shell sed -n 16p settings/parameters.txt | tr -d -c 0-9.)
VISC=$(shell sed -n 17p settings/parameters.txt | tr -d -c 0-9.-)
EPSILON_ETA=$(shell sed -n 18p settings/parameters.txt | tr -d -c 0-9.)
EPSILON_R=$(shell sed -n 19p settings/parameters.txt | tr -d -c 0-9.)

#Compiling parameters
CXX=mpic++
FLAGS=-std=c++11 -O3 $(MFEM_FLAGS)
RUN=mpirun -np $(PROCCESORS) ./
SOURCES=$(wildcard code/*.cpp)
DEPENDENCIES=$(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all mesh graph clean oclean

all: main.x results/mesh.msh
	@echo -e 'Running program ... \n'
	@$(RUN)$< --mesh results/mesh.msh --height $(HEIGHT) --internal-radius $(INNER_RADIUS) --outer-radius $(OUTER_RADIUS) \
  					--order $(ORDER) --refinements $(REF) \
					--temperature_fusion $(T_FU) --DeltaT $(DELTA_T) \
					--viscosity $(VISC) --epsilon_eta $(EPSILON_ETA) --epsilon_r $(EPSILON_R)
	@echo -e '\n'
	@cat results/convergence.txt
	@echo -e '\nDone!'

mesh: results/mesh.msh
	@echo 'Mesh created.'

graph:
ifeq ($(SHARE_DIR), NULL)
	@echo 'No share directory.'
else 
	@echo -e 'Moving graphs ... \c'
	@rm -rf $(SHARE_DIR)/step-3
	@cp -r results/graph $(SHARE_DIR)/step-3
	@echo 'Done!'
endif

main.x: $(DEPENDENCIES)
	@echo -e 'Compiling' $@ '... \c'
	@$(CXX) $(FLAGS) $^ $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

.objects/%.o: code/%.cpp
	@echo -e 'Building' $@ '... \c'
	@$(CXX) $(FLAGS) -c $< $(MFEM_LIBS) -o $@
	@echo -e 'Done!\n'

results/mesh.msh: settings/parameters.txt
	@echo -e 'Reading parameters ... \c'
	@bash settings/configure_script.sh
	@echo -e 'Generating mesh ... \c'
	@gmsh $(SCRIPT) -format msh2 -o results/mesh.msh -3 > /dev/null
	@echo -e 'Done!\n'

clean:
	@rm -rf *.x results/graph/* results/*.txt

oclean:
	@rm -rf .objects/*.o
