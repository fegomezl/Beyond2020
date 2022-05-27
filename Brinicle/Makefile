include ../local_config.mk

#Mesh parmeters
SCRIPT = $(shell sed -n 11p settings/parameters.txt | cut -d '#' -f 1)
R      = $(shell sed -n 12p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
Z      = $(shell sed -n 13p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
RIn    = $(shell sed -n 14p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)

#Simulation parameters
DT   = $(shell sed -n 17p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
T_FI = $(shell sed -n 18p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
VIS  = $(shell sed -n 19p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)

#FE parameters
REF     = $(shell sed -n 22p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
ORDER   = $(shell sed -n 23p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
ABST_C  = $(shell sed -n 24p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
RELT_C  = $(shell sed -n 25p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
ITER_C  = $(shell sed -n 26p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
ABST_S  = $(shell sed -n 27p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
RELT_S  = $(shell sed -n 28p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
EPSILON = $(shell sed -n 29p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)

#Brinicle conditions
Q  = $(shell sed -n 32p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
Nl = $(shell sed -n 33p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
Nh = $(shell sed -n 34p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
Ti = $(shell sed -n 35p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.-)
To = $(shell sed -n 36p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.-)
Tn = $(shell sed -n 37p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.-)
Si = $(shell sed -n 38p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
So = $(shell sed -n 39p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)
Sn = $(shell sed -n 40p settings/parameters.txt | cut -d '#' -f 1 | tr -d -c 0-9.)

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
	@$(RUN)$< --mesh results/mesh.msh \
			  -R $(R) -Z $(Z) -R_in $(RIn) \
			  -dt $(DT) -t_f $(T_FI) -v_s $(VIS) \
			  -ref $(REF) -o $(ORDER) \
			  -abstol_c $(ABST_C) -reltol_c $(RELT_C) -iter_c $(ITER_C) \
			  -abstol_s $(ABST_S) -reltol_s $(RELT_S) -eps $(EPSILON) \
			  -q $(Q) -nl $(Nl) -nh $(Nh) \
			  -Ti $(Ti) -To $(To) -Tn $(Tn) \
			  -Si $(Si) -So $(So) -Sn $(Sn) 
	@echo -e '\nDone!\n'

mesh: results/mesh.msh
	@echo 'Mesh created.'

graph:
ifeq ($(SHARE_DIR), NULL)
	@echo 'No share directory.'
else
	@echo -e 'Moving graphs ... \c'
	@rm -rf $(SHARE_DIR)/brinicle
	@cp -r results/graph $(SHARE_DIR)/brinicle
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

results/mesh.msh: settings/parameters.txt settings/mesh_config.sh
	@echo -e 'Reading parameters ... \c'
	@bash settings/mesh_config.sh
	@echo -e 'Generating mesh ... \c'
	@gmsh $(SCRIPT) -format msh2 -o results/mesh.msh -3 > /dev/null
	@echo -e 'Done!\n'

clean:
	@rm -rf *.x results/graph/*

oclean:
	@rm -rf .objects/*.o
