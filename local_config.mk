#Computer parameters
MFEM_INSTALL_DIR = ~/Repos/spack/opt/spack/linux-arch-sandybridge/gcc-11.1.0/mfem-develop-ih7axzrrz5if4aau5gfee62i2u6o5ham

SHARE_DIR = /media/sf_Archdata/
PROCCESORS = 4

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
