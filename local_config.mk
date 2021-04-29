#Computer parameters
MFEM_INSTALL_DIR = /homes/orionsan/cdelv/MFEM_install/mfem-4.2/config/config.mk
SHARE_DIR = /homes/orionsan/cdelv/brinicle/Share
PROCCESORS = 8

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR) #/share/mfem/config.mk
include $(CONFIG_MK)
