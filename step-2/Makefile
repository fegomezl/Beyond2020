#Computer parameters
MFEM_INSTALL_DIR = ~/spackgcc/gcc-10.2.0/mfem-4.2.0-vxsy2qfd4okzgqnuybubclbv2fnntoga/
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
SHARE_DIR = /media/sf_Shared_Folder_Arch2/
PROCCESORS = 4

include $(CONFIG_MK)

#Mesh parameters
NMESH = $(shell sed -n 2p data/parameters.txt)  
MESH = $(filter-out #mesh file,$(NMESH))
OUTER_RADIUS = $(shell sed -n 4p data/parameters.txt | tr -d -c 0-9.)
INNER_RADIUS = $(shell sed -n 5p data/parameters.txt | tr -d -c 0-9.)

#Program parameters
REF = $(shell sed -n 12p data/parameters.txt | tr -d -c 0-9.) 
ORDER = $(shell sed -n 13p data/parameters.txt | tr -d -c 0-9.)

CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS)
RUN = mpirun -np $(PROCCESORS) ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all graph clean oclean

all: main.x data/mesh_script.geo data/parameters.txt
	@echo -e 'Reading parameters ... \c'
	@bash script/parameters_script.sh
	@echo -e 'Generating mesh ... \n'
	@gmsh data/mesh_script.geo -format msh2 -o $(MESH) -3
	@echo -e '\nDone!\n'
	@echo -e 'Running program ... \n'
	@$(RUN)$< --mesh $(MESH) --order $(ORDER) --refinements $(REF) 
	@echo -e '\nDone!\n'

graph:
	@echo -e 'Moving graphs ... \c'
	@cp -r graph/ $(SHARE_DIR)
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
