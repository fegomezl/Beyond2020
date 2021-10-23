include ../../local_config.mk

#Mesh parmeters
SCRIPT=$(shell sed -n 2p settings/parameters.txt | cut -d '#' -f 1)
RMIN=$(shell sed -n 3p settings/parameters.txt | tr -d -c 0-9.-)
RMAX=$(shell sed -n 4p settings/parameters.txt | tr -d -c 0-9.-)
ZMIN=$(shell sed -n 5p settings/parameters.txt | tr -d -c 0-9.-)
ZMAX=$(shell sed -n 6p settings/parameters.txt | tr -d -c 0-9.-)

#Simulation parameters
DT=$(shell sed -n 11p settings/parameters.txt | tr -d -c 0-9.)
T_FI=$(shell sed -n 12p settings/parameters.txt | tr -d -c 0-9.)
VIS=$(shell sed -n 13p settings/parameters.txt | tr -d -c 0-9.)

#FE parameters
REF=$(shell sed -n 16p settings/parameters.txt | tr -d -c 0-9.)
ORDER=$(shell sed -n 17p settings/parameters.txt | tr -d -c 0-9.)
ABST_C=$(shell sed -n 18p settings/parameters.txt | tr -d -c 0-9.)
RELT_C=$(shell sed -n 19p settings/parameters.txt | tr -d -c 0-9.)
ITER_C=$(shell sed -n 20p settings/parameters.txt | tr -d -c 0-9.)
ABST_S=$(shell sed -n 21p settings/parameters.txt | tr -d -c 0-9.)
RELT_S=$(shell sed -n 22p settings/parameters.txt | tr -d -c 0-9.)

#Physical parameters
A_L=$(shell sed -n 25p settings/parameters.txt | tr -d -c 0-9.)
A_S=$(shell sed -n 26p settings/parameters.txt | tr -d -c 0-9.)
D_L=$(shell sed -n 27p settings/parameters.txt | tr -d -c 0-9.)
D_S=$(shell sed -n 28p settings/parameters.txt | tr -d -c 0-9.)
L_L=$(shell sed -n 29p settings/parameters.txt | tr -d -c 0-9.)
L_S=$(shell sed -n 30p settings/parameters.txt | tr -d -c 0-9.)
DELTA_T=$(shell sed -n 31p settings/parameters.txt | tr -d -c 0-9.)
EPSILON_T=$(shell sed -n 32p settings/parameters.txt | tr -d -c 0-9.)

#Compiling parameters
CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS)
RUN = mpirun -np $(PROCCESORS) ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all main mesh graph clean oclean

all: results/mesh.msh main

main: main.x
	@echo -e 'Running program ... \n'
	@$(RUN)$< --mesh results/mesh.msh -Rmin $(RMIN) -Rmax $(RMAX) -Zmin $(ZMIN) -Zmax $(ZMAX) \
			  -dt $(DT) -t_f $(T_FI) -v_s $(VIS) \
			  -ref $(REF) -o $(ORDER) \
			  -abstol_c $(ABST_C) -reltol_c $(RELT_C) -iter_c $(ITER_C) \
			  -abstol_s $(ABST_S) -reltol_s $(RELT_S) \
			  -a_l $(A_L) -a_s $(A_S) \
			  -d_l $(D_L) -d_s $(D_S) \
			  -L_l $(L_L) -L_s $(L_S) \
			  -DT $(DELTA_T) -ET $(EPSILON_T)
	@echo -e '\nDone!\n'

mesh: results/mesh.msh
	@echo 'Mesh created.'

graph:
ifeq ($(SHARE_DIR), NULL)
	@echo 'No share directory.'
else
	@echo -e 'Moving graphs ... \c'
	@rm -rf $(SHARE_DIR)/double_diffusion
	@cp -r results/graph $(SHARE_DIR)/double_diffusion
	@echo -e 'Done!'
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
	@rm -rf *.x results/graph/*

oclean:
	@rm -rf .objects/*.o
