#Computer parameters
MFEM_INSTALL_DIR = ~/spack/opt/spack/linux-debian9-skylake_avx512/gcc-11.1.0/mfem-4.2.0-ddfslg7jyfab423odotq2qqeztywajbf

SHARE_DIR = 
PROCCESORS = 32

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
