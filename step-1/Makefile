#computer parameters
MFEM_INSTALL_DIR = ~/spackgcc/gcc-10.2.0/mfem-4.2.0-vxsy2qfd4okzgqnuybubclbv2fnntoga/
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
SHARE_DIR = /media/sf_Shared_Folder_Arch2/
PROCESSORS = 4

include $(CONFIG_MK)

#Mesh parameters
MESH = data/transfinite_cylinder.msh

#Program parameters
REF = 3
ORDER = 1

CXX = mpic++
FLAGS = -std=c++11 -g -O3 $(MFEM_FLAGS)
RUN = mpirun -np 4 ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all graph clean oclean

all: main.x data/mesh_script.geo
	@echo -e 'Generating mesh ... \n'
	@gmsh data/mesh_script.geo -format msh2 -o $(MESH) -3
	@echo -e '\nDone!\n'
	@echo -e 'Running program ... \n'
	@$(RUN)$< -o $(ORDER) -r $(REF) -m $(MESH)
	@echo -e '\n'
	@cat data/convergence.txt
	@echo -e '\nDone!\n'
graph:
	@echo -e 'Moving graphs ... \c'
	@cp -r graph/ $(SHARE_DIR)
	@echo -e 'Done\n'

main.x: $(DEPENDENCIES)
	@echo -e 'Compiling' $@ '... \c'
	@$(CXX) $(FLAGS) $^ $(MFEM_LIBS) -o $@
	@echo -e 'Done\n'

.objects/%.o: code/%.cpp
	@echo -e 'Building' $@ '... \c'
	@$(CXX) $(FLAGS) -c $< $(MFEM_LIBS) -o $@
	@echo -e 'Done\n'

clean:
	@rm -rf *.x graph/*

oclean:
	@rm -rf .objects/*.o
