include ../local_config.mk
include settings/parameters.mk

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

results/mesh.msh: settings/parameters.mk
	@echo -e 'Generating mesh ... \c'
	@$(GMSH_INSTALL)gmsh $(SCRIPT) -setnumber R $$(( $(R)/$(RIn) )) -setnumber Z $$(( $(Z)/$(RIn) )) -setnumber r 1 -format msh2 -o results/mesh.msh -3 > /dev/null
	@echo -e 'Done!\n'

clean:
	@rm -rf *.x results/graph/*

oclean:
	@rm -rf .objects/*.o
