MFEM_INSTALL_DIR = ~/spackgcc/gcc-10.2.0/mfem-4.2.0-vxsy2qfd4okzgqnuybubclbv2fnntoga/
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
MESH = data/transfinite_cylinder.msh
REF = 0
ORDER = 1

include $(CONFIG_MK)

CXX = mpic++
FLAGS = -std=c++11 -g -O3 $(MFEM_FLAGS)
RUN = mpirun -np 4 ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all graph clean oclean

all: main.x
	@echo -e 'Running program ... \n'
	@$(RUN)$< -o $(ORDER) -r $(REF) -m $(MESH)
	@echo -e '\nDone\n'

graph:
	@echo -e 'Moving graphs ... \c'
	@cp -r graph/ /media/sf_Shared_Folder_Arch2/
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
