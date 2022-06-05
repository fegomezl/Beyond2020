# --------------------------
# Parameters of the program
# --------------------------
# Units:
# mm      (milimeters)
# s       (seconds)
# °C      (celcius)
# %wt     (weight_percentage)

#Mesh Parameters:
SCRIPT = settings/sea.geo    #Mesh_script
R = 10                  #R(mm)
Z = 40                  #Z(mm)
RIn = 2                   #Inflow(mm)

#Simulation parameters:
DT = 0.001         #Dt(s)
T_FI = 1           #Final_time(s)
VIS = 5            #Visualization_steps

#FE parameters:
REF = 4          #Refinements
ORDER = 1          #Order
ABST_C= 0          #abstol(Conduction)
RELT_C = 0.00000001 #reltol(Conduction)
ITER_C = 100        #iter(Conduction)
ABST_S = 0.00001    #abstol(SUNDIALS)
RELT_S = 0.00001    #reltol(SUNDIALS)
EPSILON = 9          #nEpsilon

#Brinicle conditions:
Q = 300        #FlowRate(mm^3/s)
Nl = 1          #Nucleation_length(mm)
Nh = 4          #Nucleation_height(mm)
Ti = -2         #Initial_temperature(°C)
To = -15         #Inflow_temperature(°C)
Tn = -10         #Nucleation_temperature(°C)
Si = 3.5        #Initial_salinity(%wt)
So = 22.5       #Inflow_salinity(%wt)
Sn = 3.5        #Nucleation_salinity(%wt)
