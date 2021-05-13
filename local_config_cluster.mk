#Computer parameters
MFEM_INSTALL_DIR = ~/MFEM_install/mfem-4.2/
SHARE_DIR = 
PROCCESORS = 8

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/config/config.mk
include $(CONFIG_MK)
