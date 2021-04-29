#Computer parameters
MFEM_INSTALL_DIR =/home/wind/Library/spack/opt/spack/linux-pop20-excavator/gcc-10.2.0/mfem-4.2.0-ie7sv2es2hutigzmv6ppypyrot2vx7k5
SHARE_DIR =/home/wind/Pictures/Beyond
PROCCESORS = 4

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
