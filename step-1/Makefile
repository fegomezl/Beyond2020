MFEM_INSTALL_DIR = ~/repos/spack/opt/spack/linux-arch-sandybridge/gcc-10.2.0/mfem-4.2.0-cnxxrnnnbzf2o6blvzgu3g4546jmf2ku
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
MESH = data/cylinder_hollow.msh

include $(CONFIG_MK)

CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS)
RUN = mpirun -np 4 ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all graph clean oclean

all: main.x
	@echo -e 'Running program ... \n'
	@$(RUN)$< -o 1 -r 4 
	@echo -e '\nDone\n'

graph:
	@echo -e 'Moving graphs ... \c'
	@cp -r graph/ /media/sf_Arch_data/
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
