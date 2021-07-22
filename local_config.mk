#Computer parameters
MFEM_INSTALL_DIR = /opt/spack/opt/spack/linux-pop21-icelake/gcc-10.3.0/mfem-develop-5oquzzvtgmiaq3iwn57wj2opviva6wxd


SHARE_DIR = 
PROCCESORS = 8

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
