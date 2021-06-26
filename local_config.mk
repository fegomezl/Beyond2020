#Computer parameters
#MFEM_INSTALL_DIR = ~/Library/spack/opt/spack/linux-pop20-excavator/gcc-10.2.0/mfem-4.2.0-linjgmwmtphlcfh4ow6hopjpq4aw47oz

MFEM_INSTALL_DIR =/home/wind/Library/MFEM/mfem/build/config
SHARE_DIR = 
PROCCESORS = 4

#Add variables from MFEM
#CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
CONFIG_MK = $(MFEM_INSTALL_DIR)/config.mk
include $(CONFIG_MK)
