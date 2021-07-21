#Computer parameters
MFEM_INSTALL_DIR = ~/spackgcc/mfem-4.2.0-memory-superlu-ck7pgi4afqck7k6dn2w4cisxfk7db3l6/
SHARE_DIR = /media/sf_Shared_Folder_Arch2/
PROCCESORS = 4

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
