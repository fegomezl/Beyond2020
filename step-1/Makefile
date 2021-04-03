MFEM_INSTALL_DIR = ~/spackgcc/gcc-10.2.0/mfem-4.2.0-vxsy2qfd4okzgqnuybubclbv2fnntoga/
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
MESH = data/star.mesh

include $(CONFIG_MK)

CXX = mpic++
FLAGS = -std=c++11 -O3 $(MFEM_FLAGS)
RUN = mpirun -np 4 ./
SOURCES = $(wildcard code/*.cpp)
DEPENDENCIES = $(SOURCES:code/%.cpp=.objects/%.o)

.PHONY: all graph clean oclean

all: main.x
	@echo -e 'Running program ... \n'
<<<<<<< HEAD
	@$(RUN)$< -m $(MESH)
=======
	@$(RUN)$< -o 1 -r 3 
>>>>>>> f4c7f9434a7405927babf912e0268e10b27c4a3c
	@echo -e '\nDone\n'

graph:
	@echo -e 'Moving graphs ... \c'
	@cp -rf graph/ /media/sf_Shared_Folder_Arch2/
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
