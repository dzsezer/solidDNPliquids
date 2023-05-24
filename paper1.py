####
# Author: Deniz Sezer
# Paper I: Dynamic view of the solid-state DNP effect
# https://doi.org/10.5194/mr-2023-1
###
# Code for Figure 5, 6b, 7 8, 9 and 10
###

import plots

# Microwave magnetic fields in Gauss
B1s = [0.5,1.5,3,6] # G


# Electronic relaxation times in micro seconds
T2e = 60e-3 #us
T1e = 9*T2e 

# Figure 5 
plots.plot_bloch(B1s,T1e,T2e)


freqIs = [14,25,45,140]
sides  = [16,20,30,50]

# Figure 6b
plots.plot_Fxyz(B1s,freqIs,sides,T1e,T2e)


#Microwave magnetic field in Gauss
B1 = 6 #G

# Figure 7
plots.plot_newbloch(B1,freqIs,sides,T1e,T2e)
# Figure 8
plots.plot_rates(B1,freqIs,sides,T1e,T2e)


# Nuclear T1 in micro seconds
T1n = 30e3 # us
#distance of closest approach
closest_distance_nm = 1.0 #nm
#radical concentration
radical_conc_molar = 0.1 #molar

#Figure 9
plots.plot_enhancement(B1,freqIs,sides,T1e,T2e,T1n,closest_distance_nm,radical_conc_molar)

#Figure 10a
plots.plot_filters(B1,freqIs,sides,T1e,T2e,'class')
#Figure 10b
plots.plot_filters(B1,freqIs,sides,T1e,T2e,'new')

