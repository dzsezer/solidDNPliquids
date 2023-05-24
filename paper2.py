####
# Author: Deniz Sezer
# Paper II: The solid-state DNP effect in liquids
# https://doi.org/10.5194/mr-2023-2
###
# Code for Figure 3, 4, 5, 6, 7 and 8
###
import plots, funs

B1s = [0.5,1.5,3,6] # G

T2e = 60e-3 #micro s
T1e = 9*T2e

tau = T2e/5


T1n_min = 50e3 # micro s
#distance of closest approach
closest_distance_nm = 1. #nm
#radical concentration
radical_conc_molar = 0.1 #molar

#R1n = 1/T1n
delta2 = funs.delta2(closest_distance_nm, radical_conc_molar)

freqIs = [14,45,140,400]
sides  = [32,60,100,140]


B1 = 6 #G

# Figure 3
plots.plot_liquid('solid',B1,freqIs,sides,T1e,T2e,tau,delta2,T1n_min)
# Figure 4
plots.plot_liquid('liquid',B1,freqIs,sides,T1e,T2e,tau,delta2,T1n_min)
# Figure 5
plots.plot_liquid('pert',B1,freqIs,sides,T1e,T2e,tau,delta2,T1n_min)


# Figure 6
plots.plot_Tx(freqIs[-1],sides[-1],B1,T1e,T2e,tau)

#############################
# Comparison with experiment
#####################


#read experimental EPR spectrum
with open('data/figureS8-EPR.dat') as f:
    lines = f.readlines()
epr_J = funs.process_2columns(lines[1:])

#manual fit to EPR spectrum to determine electronic T2
T2e = 215e-3 # micro s
phase = -9 # degrees
plots.plot_epr(epr_J,T2e,phase,scale=98)

#read experimental DNP profile
with open('data/figureS8-DNP.dat') as f:
    lines = f.readlines()
dnp_J = funs.process_2columns(lines[1:])

########################
# Set B1 and nuclear T1 
########################
T1n = 50e3 # micro s
R1n = 1/T1n

B1 = 6 # Gauss

#######################################
# Specify initial values for the FIT
########################################
closest_distance_nm = 1.0 #nm
radical_conc_molar = 0.1 #molar

T1e = 8.6*T2e
tau = T2e/30.4

fit = True

# Figure 7
plots.plot_dnp(dnp_J,B1,T2e,T1e,R1n,
	closest_distance_nm,radical_conc_molar,tau,'ffhs',fit)

# Figure 8
plots.plot_dnp(dnp_J,B1,T2e,T1e,R1n,
	closest_distance_nm,radical_conc_molar,tau,'exp',fit)




