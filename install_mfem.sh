#Exit when any command fails
set -e
#Print commands 
set -x   

INSTALL_DIR=$(pwd)
PETSC_DIR=${INSTALL_DIR}/petsc
SUNDIALS_DIR=${INSTALL_DIR}/sundials
MFEM_INSTALL_DIR=${INSTALL_DIR}/mfem/build

####################
### Download Source
####################

#Get Sundials
wget https://github.com/LLNL/sundials/releases/download/v5.8.0/sundials-5.8.0.tar.gz
tar -xf sundials-5.8.0.tar.gz
mv sundials-5.8.0 sundials

#Get Petsc
git clone -b release https://gitlab.com/petsc/petsc.git petsc_src

#Get MFEM
git clone https://github.com/mfem/mfem.git

####################
### INSTALL PETCs
####################

cd ${INSTALL_DIR}/petsc_src

./configure --prefix=${INSTALL_DIR}/petsc --with-debugging=0 --download-hypre --download-superlu_dist --download-metis --download-parmetis --download-slepc --download-sprng --with-fortran-bindings=0 COPTFLAGS="-O3 -march=native -mtune=native" CXXOPTFLAGS="-O3 -march=native -mtune=native" FOPTFLAGS="-O3 -march=native -mtune=native"

#In case you dont have MPI on your system add the following line to above command
#--download-mpich  --download-hwloc
#
#In case you dont have zlib on your system add the following line to above command
#--download-zlib
#
#If you have MPI or zlib adding these lines will cause truble
#
#If you are runing on a virtual variable/module base system add after ./configure
# CC=$CC CXX=$CXX FC=$FC F77=$F77 F90=$F90 CPP=$CPP 

make PETSC_DIR=${INSTALL_DIR}/petsc_src PETSC_ARCH=arch-linux-c-opt all
make PETSC_DIR=${INSTALL_DIR}/petsc_src PETSC_ARCH=arch-linux-c-opt install
make SLEPC_DIR=${INSTALL_DIR}/petsc PETSC_DIR=${INSTALL_DIR}/petsc PETSC_ARCH="" check

#If you want the the system to use the MPI executable and/or the installed librarys as default add them to the path
#export PATH=${INSTALL_DIR}/petsc/bin/:$PATH
#export LD_LIBRARY_PATH=${INSTALL_DIR}/petsc/lib:$LD_LIBRARY_PATH

####################
### INSTALL SUNDAILS
####################

cd ${SUNDIALS_DIR}
mkdir build
cd ${SUNDIALS_DIR}/build
cmake -DCMAKE_INSTALL_PREFIX=${SUNDIALS_DIR}/install -DEXAMPLES_INSTALL_PATH=${SUNDIALS_DIR}/install/examples -DENABLE_MPI:BOOL=ON -DENABLE_PETSC:BOOL=ON -DENABLE_HYPRE:BOOL=ON -DHYPRE_INCLUDE_DIR=${PETSC_DIR}/include -DHYPRE_LIBRARY_DIR=${PETSC_DIR}/lib DLAPACK_LIBRARIES=${PETSC_DIR}/lib -DPETSC_ARCH= -DPETSC_DIR=${PETSC_DIR} -DSUNDIALS_INDEX_SIZE=32 ${SUNDIALS_DIR}

make install -j 8

#If RAM problems when executing make, add 4 or 8 afte the -j

####################
### INSTALL MFEM
####################

cd ${INSTALL_DIR}/mfem
mkdir build
cd ${INSTALL_DIR}/mfem/build
cmake -DMFEM_USE_MPI:BOOL=ON -DMFEM_USE_METIS:BOOL=ON -DMFEM_ENABLE_MINIAPPS:BOOL=ON -DMFEM_USE_SUNDIALS:BOOL=ON -DMFEM_USE_SLEPC:BOOL=ON -DMFEM_USE_PETSC:BOOL=ON -DMFEM_USE_ZLIB:BOOL=ON -DMFEM_USE_SUPERLU5:BOOL=ON -DHYPRE_DIR=${PETSC_DIR} -DMETIS_DIR=${PETSC_DIR} -DParMETIS_DIR=${PETSC_DIR} -DPETSC_ARCH= -DPETSC_DIR=${PETSC_DIR} -DPETSC_EXECUTABLE_RUNS=${PETSC_DIR}/lib/petsc/bin -DSUNDIALS_DIR=${SUNDIALS_DIR}/install -DSLEPC_DIR=${PETSC_DIR} -DSLEPC_VERSION_OK=yes -DMFEM_USE_SUPERLU=yes -DSuperLUDist_DIR=${PETSC_DIR} ${INSTALL_DIR}/mfem -DCMAKE_INSTALL_NAME_DIR=${MFEM_INSTALL_DIR}/lib -DCMAKE_INSTALL_PREFIX=${MFEM_INSTALL_DIR} -DCMAKE_INSTALL_RPATH=${MFEM_INSTALL_DIR}/lib

#In case you are using the zlib instaled with PETSC add before ${INSTALL_DIR}/mfem
#-DZLIB_INCLUDE_DIR=${PETSC_DIR}/include -DZLIB_LIBRARIES=${PETSC_DIR}/lib/libz.a

make -j 8
make examples -j 8
make miniapps -j 8
make tests -j 8
make install -j 8
make test

#If RAM problems when executing make, add 4 or 8 afte the -j