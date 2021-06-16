#Computer parameters
MFEM_INSTALL_DIR = ~/Repos/spack/opt/spack/linux-arch-sandybridge/gcc-11.1.0/mfem-4.2.0-4ybpcz34b6xk64a3acnu4u5lnkd3a5jr

SHARE_DIR = /media/sf_Archdata/
PROCCESORS = 4

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
