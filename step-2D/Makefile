include ../local_config.mk

CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS)
RUN = mpirun -np $(PROCCESORS) ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)
MESH = graph/mesh.msh
SCRIPT = scripts/structured_square.geo

.PHONY: all graph clean oclean

all: main.x
	@echo -e 'Generating mesh ... \n'
	@gmsh $(SCRIPT) -format msh2 -o $(MESH) -3 > /dev/null
	@echo -e '\nDone!\n'
	@echo -e 'Running program ... \n'
	@$(RUN)$< --mesh $(MESH)
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
