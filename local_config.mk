#Computer parameters
MFEM_INSTALL_DIR = ~/spackgcc/mymfem-4.2.0-h4kgyir56to4diecuv6jzvwzpurmxbec/
SHARE_DIR = /media/sf_Shared_Folder_Arch2/
PROCCESORS = 4

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
