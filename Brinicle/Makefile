include ../local_config.mk

#Mesh parmeters
SCRIPT=$(shell sed -n 2p settings/parameters.txt | cut -d '#' -f 1)
RMIN=$(shell sed -n 3p settings/parameters.txt | tr -d -c 0-9.-)
RMAX=$(shell sed -n 4p settings/parameters.txt | tr -d -c 0-9.-)
ZMIN=$(shell sed -n 5p settings/parameters.txt | tr -d -c 0-9.-)
ZMAX=$(shell sed -n 6p settings/parameters.txt | tr -d -c 0-9.-)
LI=$(shell sed -n 7p settings/parameters.txt | tr -d -c 0-9.-)
LO=$(shell sed -n 8p settings/parameters.txt | tr -d -c 0-9.-)

#Simulation parameters
DT=$(shell sed -n 13p settings/parameters.txt | tr -d -c 0-9.)
T_FI=$(shell sed -n 14p settings/parameters.txt | tr -d -c 0-9.)
VIS=$(shell sed -n 15p settings/parameters.txt | tr -d -c 0-9.)
RC=$(shell sed -n 16p settings/parameters.txt | tr -d -c 0-9)

#FE parameters
REF=$(shell sed -n 19p settings/parameters.txt | tr -d -c 0-9.)
ORDER=$(shell sed -n 20p settings/parameters.txt | tr -d -c 0-9.)
ABST_C=$(shell sed -n 21p settings/parameters.txt | tr -d -c 0-9.)
RELT_C=$(shell sed -n 22p settings/parameters.txt | tr -d -c 0-9.)
ITER_C=$(shell sed -n 23p settings/parameters.txt | tr -d -c 0-9.)
ABST_S=$(shell sed -n 24p settings/parameters.txt | tr -d -c 0-9.)
RELT_S=$(shell sed -n 25p settings/parameters.txt | tr -d -c 0-9.)

#Physical parameters
T_FU=$(shell sed -n 28p settings/parameters.txt | tr -d -c 0-9.-)
DELTA_T=$(shell sed -n 29p settings/parameters.txt | tr -d -c 0-9.)
EPSILON_T=$(shell sed -n 30p settings/parameters.txt | tr -d -c 0-9.)
EPSILON_ETA=$(shell sed -n 31p settings/parameters.txt | tr -d -c 0-9.)
C_L=$(shell sed -n 32p settings/parameters.txt | tr -d -c 0-9.)
C_S=$(shell sed -n 33p settings/parameters.txt | tr -d -c 0-9.)
K_L=$(shell sed -n 34p settings/parameters.txt | tr -d -c 0-9.)
K_S=$(shell sed -n 35p settings/parameters.txt | tr -d -c 0-9.)
D_L=$(shell sed -n 36p settings/parameters.txt | tr -d -c 0-9.)
D_S=$(shell sed -n 37p settings/parameters.txt | tr -d -c 0-9.)
L_L=$(shell sed -n 38p settings/parameters.txt | tr -d -c 0-9.)
L_S=$(shell sed -n 39p settings/parameters.txt | tr -d -c 0-9.)

#Brinicle conditions
V=$(shell sed -n 42p settings/parameters.txt | tr -d -c 0-9.)
Ti=$(shell sed -n 43p settings/parameters.txt | tr -d -c 0-9.-)
To=$(shell sed -n 44p settings/parameters.txt | tr -d -c 0-9.-)
Si=$(shell sed -n 45p settings/parameters.txt | tr -d -c 0-9.)
So=$(shell sed -n 46p settings/parameters.txt | tr -d -c 0-9.)
Nl=$(shell sed -n 47p settings/parameters.txt | tr -d -c 0-9.)
Nh=$(shell sed -n 48p settings/parameters.txt | tr -d -c 0-9.)
Tn=$(shell sed -n 49p settings/parameters.txt | tr -d -c 0-9.-)
Sn=$(shell sed -n 50p settings/parameters.txt | tr -d -c 0-9.)

#Restart conditions
R=$(shell sed -n 53p settings/parameters.txt | tr -d -c 0-9.)
T_IN=$(shell sed -n 54p settings/parameters.txt | tr -d -c 0-9.)

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
						-Li $(LI) -Lo $(LO) \
						-dt $(DT) -t_f $(T_FI) -v_s $(VIS) \
						-ref $(REF) -o $(ORDER) \
						-abstol_c $(ABST_C) -reltol_c $(RELT_C) -iter_c $(ITER_C) \
						-abstol_s $(ABST_S) -reltol_s $(RELT_S) \
						-T_f $(T_FU) -DT $(DELTA_T) -ET $(EPSILON_T) -EEta $(EPSILON_ETA) \
						-c_l $(C_L) -c_s $(C_S) -k_l $(K_L) -k_s $(K_S) \
						-d_l $(D_L) -d_s $(D_S) -L_l $(L_L) -L_s $(L_S) \
						-v $(V) -Ti $(Ti) -To $(To) -Si $(Si) -So $(So) \
						-nl $(Nl) -nh $(Nh) -Tn $(Tn) -Sn $(Sn) \
						-r $(R) -t_i $(T_IN)
	@sed -i 54d settings/parameters.txt
	@echo -e '\nDone!\n'

mesh: results/mesh.msh
	@echo 'Mesh created.'

dt: results/analysis.pdf
	@xpdf $<

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

results/mesh.msh: settings/parameters.txt
	@echo -e 'Reading parameters ... \c'
	@bash settings/configure_script.sh
	@echo -e 'Generating mesh ... \c'
	@gmsh $(SCRIPT) -format msh2 -o results/mesh.msh -3 > /dev/null
	@echo -e 'Done!\n'

results/analysis.pdf: settings/analysis.gp
	@gnuplot settings/analysis.gp

clean:
	@rm -rf *.x results/graph/*

rclean:
	@rm -rf *.x results/restart/*.msh results/restart/*.gf results/*.txt results/*.pdf
	@echo '0          #Initial_time' >> settings/parameters.txt
	@sed -i 54d settings/parameters.txt

oclean:
	@rm -rf .objects/*.o
