include ../../local_config.mk

#Mesh parmeters
SCRIPT=$(shell sed -n 2p settings/parameters.txt | cut -d '#' -f 1)
RMIN=$(shell sed -n 3p settings/parameters.txt | tr -d -c 0-9.-)
RMAX=$(shell sed -n 4p settings/parameters.txt | tr -d -c 0-9.-)
ZMIN=$(shell sed -n 5p settings/parameters.txt | tr -d -c 0-9.-)
ZMAX=$(shell sed -n 6p settings/parameters.txt | tr -d -c 0-9.-)

#Program parameters
REF_I = $(shell sed -n 11p settings/parameters.txt | tr -d -c 0-9.)
REF_L = $(shell sed -n 12p settings/parameters.txt | tr -d -c 0-9.)
ORDER = $(shell sed -n 13p settings/parameters.txt | tr -d -c 0-9.)
T_FU = $(shell sed -n 16p settings/parameters.txt | tr -d -c 0-9.-)
C_L = $(shell sed -n 17p settings/parameters.txt | tr -d -c 0-9.)
C_S = $(shell sed -n 18p settings/parameters.txt | tr -d -c 0-9.)
K_L = $(shell sed -n 19p settings/parameters.txt | tr -d -c 0-9.)
K_S = $(shell sed -n 20p settings/parameters.txt | tr -d -c 0-9.)
L = $(shell sed -n 21p settings/parameters.txt | tr -d -c 0-9.)
DELTA_T_I = $(shell sed -n 22p settings/parameters.txt | tr -d -c 0-9.)
DELTA_T_JS = $(shell sed -n 23p settings/parameters.txt | tr -d -c 0-9.)
DELTA_T_J = $(shell sed -n 24p settings/parameters.txt | tr -d -c 0-9.)
DT = $(shell sed -n 25p settings/parameters.txt | tr -d -c 0-9.)
T_FI = $(shell sed -n 26p settings/parameters.txt | tr -d -c 0-9.)
VIS = $(shell sed -n 27p settings/parameters.txt | tr -d -c 0-9.)
ODE = $(shell sed -n 28p settings/parameters.txt | tr -d -c 0-9.)
RELT = $(shell sed -n 29p settings/parameters.txt | tr -d -c 0-9.)
ABST = $(shell sed -n 30p settings/parameters.txt | tr -d -c 0-9.)
STATE = $(shell sed -n 31p settings/parameters.txt | tr -d -c 0-9.)

#Compiling parameters
CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS)
RUN = mpirun -np $(PROCCESORS) ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all main mesh send check graph cluster clean oclean

all: main.x results/mesh.msh
	@echo -e 'Running program ... \n'
	@$(RUN)$< --mesh results/mesh.msh -Rmin $(RMIN) -Rmax $(RMAX) -Zmin $(ZMIN) -Zmax $(ZMAX) \
						-o $(ORDER) -r_i $(REF_I) -r_l $(REF_L) -DT_i $(DELTA_T_I) -DT_js $(DELTA_T_JS) -DT_j $(DELTA_T_J) \
						-T_f $(T_FU) -c_s $(C_S) -c_l $(C_L) -k_s $(K_S) -k_l $(K_L) -L $(L) \
						-dt $(DT) -t_f $(T_FI) -v_s $(VIS) -ode $(ODE) -reltol $(RELT) -abstol $(ABST) -state $(STATE)
	@echo -e '\nDone!\n'

main: main.x
	@echo 'Program Compiled.'

mesh: results/mesh.msh
	@echo 'Mesh created.'

send:
	@bash settings/send.sh

check:
	@tail results/brinicle.err
	@tail results/brinicle.out

graph:
ifeq ($(SHARE_DIR), NULL)
	@echo 'No share directory.'
else
	@echo -e 'Moving graphs ... \c'
	@rm -rf $(SHARE_DIR)/twophase_benchmark
	@cp -r results $(SHARE_DIR)/twophase_benchmark
	@echo -e 'Done!'
endif

cluster:
	@bash settings/send.sh

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
	@rm -rf *.x results/*

oclean:
	@rm -rf .objects/*.o
