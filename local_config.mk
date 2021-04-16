#Computer parameters
MFEM_INSTALL_DIR = #MFEM Installation Directory 
SHARE_DIR = #Where to send the graphs output
PROCCESORS = #Number of processors

#Add variables from MFEM
CONFIG_MK = $(MFEM_INSTALL_DIR)/share/mfem/config.mk
include $(CONFIG_MK)
