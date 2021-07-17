#Computer parameters
MFEM_INSTALL_DIR = ~/spackgcc/gcc-10.2.0/mfem-4.2.0-6s3377riaaks5lfhx3zoy6nybrq7axgw/
SHARE_DIR = /media/sf_Shared_Folder_Arch2/
PROCCESORS = 4

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
