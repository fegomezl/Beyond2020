#Computer parameters
MFEM_INSTALL_DIR = ~/spack/opt/spack/linux-debian9-skylake_avx512/gcc-11.1.0/mfem-develop-yqkyh5qykk3hmel72nq5d2ytubw45uqm

SHARE_DIR = 
PROCCESORS = 32

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
