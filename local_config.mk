#Computer parameters
MFEM_INSTALL_DIR = ~/repos/spack/opt/spack/linux-arch-sandybridge/gcc-10.2.0/mfem-4.2.0-4bqzidc6xxuyfw2pfocmtt5e7eumymps
SHARE_DIR = /media/sf_Arch_data/
PROCCESORS = 4

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
