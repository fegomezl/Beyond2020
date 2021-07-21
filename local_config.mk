#Computer parameters
MFEM_INSTALL_DIR = /opt/spack/opt/spack/linux-pop21-icelake/gcc-10.3.0/mfem-develop-sxen7z2vnfukq66lkx5rx3rmluvuml2l

SHARE_DIR = 
PROCCESORS = 4

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
