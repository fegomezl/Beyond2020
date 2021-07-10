#Computer parameters
MFEM_INSTALL_DIR = /opt/spack/opt/spack/linux-pop21-icelake/gcc-10.3.0/mfem-4.2.0-j2kudzxp2l65t32h5zwqum7tjff2qqqn

SHARE_DIR = 
PROCCESORS = 4

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
