include ../../local_config.mk

#Mesh parmeters
SCRIPT=$(shell sed -n 2p settings/parameters.txt | cut -d '\#' -f 1)
RMIN=$(shell sed -n 3p settings/parameters.txt | tr -d -c 0-9.-)
RMAX=$(shell sed -n 4p settings/parameters.txt | tr -d -c 0-9.-)
ZMIN=$(shell sed -n 5p settings/parameters.txt | tr -d -c 0-9.-)
ZMAX=$(shell sed -n 6p settings/parameters.txt | tr -d -c 0-9.-)

#Simulation parameters
DT=$(shell sed -n 12p settings/parameters.txt | tr -d -c 0-9.)
T_FI=$(shell sed -n 13p settings/parameters.txt | tr -d -c 0-9.)
VIS=$(shell sed -n 14p settings/parameters.txt | tr -d -c 0-9.)

#FE parameters
REF=$(shell sed -n 17p settings/parameters.txt | tr -d -c 0-9.)
ORDER=$(shell sed -n 18p settings/parameters.txt | tr -d -c 0-9.)
ODE=$(shell sed -n 19p settings/parameters.txt | tr -d -c 0-9.)
ABST_C=$(shell sed -n 20p settings/parameters.txt | tr -d -c 0-9.)
RELT_C=$(shell sed -n 21p settings/parameters.txt | tr -d -c 0-9.)
ITER_C=$(shell sed -n 22p settings/parameters.txt | tr -d -c 0-9.)
ABST_S=$(shell sed -n 23p settings/parameters.txt | tr -d -c 0-9.)
RELT_S=$(shell sed -n 24p settings/parameters.txt | tr -d -c 0-9.)

#Physical parameters
T_FU=$(shell sed -n 27p settings/parameters.txt | tr -d -c 0-9.-)
DELTA_T=$(shell sed -n 28p settings/parameters.txt | tr -d -c 0-9.)
C_L=$(shell sed -n 29p settings/parameters.txt | tr -d -c 0-9.)
C_S=$(shell sed -n 30p settings/parameters.txt | tr -d -c 0-9.)
K_L=$(shell sed -n 31p settings/parameters.txt | tr -d -c 0-9.)
K_S=$(shell sed -n 32p settings/parameters.txt | tr -d -c 0-9.)
D_L=$(shell sed -n 33p settings/parameters.txt | tr -d -c 0-9.)
D_S=$(shell sed -n 34p settings/parameters.txt | tr -d -c 0-9.)
L=$(shell sed -n 35p settings/parameters.txt | tr -d -c 0-9.)
EPSILON_ETA=$(shell sed -n 36p settings/parameters.txt | tr -d -c 0-9.)
EPSILON_R=$(shell sed -n 37p settings/parameters.txt | tr -d -c 0-9.)

#Compiling parameters
CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS)
RUN = mpirun -np $(PROCCESORS) ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all main mesh graph clean oclean

all: main.x #results/mesh.msh
	@echo -e 'Running program ... \n'
	@$(RUN)$< --mesh results/mesh.msh -Rmin $(RMIN) -Rmax $(RMAX) -Zmin $(ZMIN) -Zmax $(ZMAX) \
						-dt $(DT) -t_f $(T_FI) -v_s $(VIS) \
						-r $(REF) -o $(ORDER) -ode $(ODE) \
						-abstol_c $(ABST_C) -reltol_c $(RELT_C) -iter_c $(ITER_C) \
						-abstol_s $(ABST_S) -reltol_s $(RELT_S) \
						-T_f $(T_FU) -DT $(DELTA_T) -c_l $(C_L) -c_s $(C_S) \
						-k_l $(K_L) -k_s $(K_S) -D_l $(D_L) -D_s $(D_S) -L $(L) \
						 -e_eta $(EPSILON_ETA) -e_r $(EPSILON_R)
	@echo -e '\nDone!\n'

main: main.x
	@echo 'Program Compiled.'

mesh: results/mesh.msh
	@echo 'Mesh created.'

graph:
ifeq ($(SHARE_DIR), NULL)
	@echo 'No share directory.'
else
	@echo -e 'Moving graphs ... \c'
	@rm -rf $(SHARE_DIR)/step-4
	@cp -r results $(SHARE_DIR)/step-4
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
