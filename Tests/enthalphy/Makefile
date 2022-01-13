include ../../local_config.mk

#Mesh parmeters
SCRIPT=$(shell sed -n 2p settings/parameters.txt | cut -d '#' -f 1)
R=$(shell sed -n 3p settings/parameters.txt | tr -d -c 0-9.)
Z=$(shell sed -n 4p settings/parameters.txt | tr -d -c 0-9.)

#Simulation parameters
DT=$(shell sed -n 9p settings/parameters.txt | tr -d -c 0-9.)
T_FI=$(shell sed -n 10p settings/parameters.txt | tr -d -c 0-9.)
VIS=$(shell sed -n 11p settings/parameters.txt | tr -d -c 0-9.)

#FE parameters
REF=$(shell sed -n 14p settings/parameters.txt | tr -d -c 0-9.)
ORDER=$(shell sed -n 15p settings/parameters.txt | tr -d -c 0-9.)
ABST_C=$(shell sed -n 16p settings/parameters.txt | tr -d -c 0-9.)
RELT_C=$(shell sed -n 17p settings/parameters.txt | tr -d -c 0-9.)
ITER_C=$(shell sed -n 18p settings/parameters.txt | tr -d -c 0-9.)
ABST_S=$(shell sed -n 19p settings/parameters.txt | tr -d -c 0-9.)
RELT_S=$(shell sed -n 20p settings/parameters.txt | tr -d -c 0-9.)

#Physical parameters
JUMPSCALE=$(shell sed -n 23p settings/parameters.txt | tr -d -c 0-9.)
A_L=$(shell sed -n 24p settings/parameters.txt | tr -d -c 0-9.)
A_S=$(shell sed -n 25p settings/parameters.txt | tr -d -c 0-9.)

#Benchmark parameters
T_L=$(shell sed -n 28p settings/parameters.txt | tr -d -c 0-9.)
T_S=$(shell sed -n 29p settings/parameters.txt | tr -d -c 0-9.)

#Compiling parameters
CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS)
RUN = mpirun --use-hwthread-cpus -np $(PROCCESORS) ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all main mesh graph plot figure clean oclean

all: main.x results/mesh.msh
	@echo -e 'Running program ... \n'
	@$(RUN)$< --mesh results/mesh.msh -R $(R) -Z $(Z) \
			  -dt $(DT) -t_f $(T_FI) -v_s $(VIS) \
			  -r $(REF) -o $(ORDER) \
			  -abstol_c $(ABST_C) -reltol_c $(RELT_C) -iter_c $(ITER_C) \
			  -abstol_s $(ABST_S) -reltol_s $(RELT_S) \
			  -nJ $(JUMPSCALE) -a_l $(A_L) -a_s $(A_S) \
			  -T_l $(T_L) -T_s $(T_S)
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
	@rm -rf $(SHARE_DIR)/enthalphy
	@cp -r results $(SHARE_DIR)/enthalphy
	@echo -e 'Done!'
endif

plot:
	@$(PARAVIEW_PATH)pvpython settings/interfase.py

figure:
	@gnuplot settings/regression.gp
	@xpdf results/regression.pdf

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
