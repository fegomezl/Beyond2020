#Computer parameters
MFEM_INSTALL_DIR =~/Library/spack/opt/spack/linux-pop20-excavator/gcc-10.2.0/mfem-4.2.0-c4o6keal7hdkfy5kdoipgjdltk3ddtda
SHARE_DIR = ~/Pictures/Beyond
PROCCESORS = 4

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
