"""
Figure 5: Astrocytes connected in networks

De Pitta model with multiple compartments, but without diffusion. 
All compartments are the same, but there is a soma. 
"""
import matplotlib      # DELETE
# matplotlib.use('agg')  # DELETE
from brian2 import *
from brian2.units import *
from brian2.units.allunits import ymole, ymol, ymolar, umolar
#import brian2.units as un
#from brian2.units.constants import * # GE

# import plot_utils as pu

set_device('cpp_standalone', directory='./codegen', build_on_run=True)  # Use fast "C++ standalone mode"

################################################################################
# Model parameters
################################################################################
### General parameters
## Something wrong with the time scales. I should be able to have much larger time steps
duration = 1000*msecond       # Total simulation time
sim_dt = 40.*ms               # Integrator/sampling step (too small!!)  (beyond 0.30*ms, there is 2-point oscillation in solution)
#duration = 50*second       # Total simulation time
#sim_dt = 1*second               # Integrator/sampling step (too small!!)  (beyond 0.30*ms, there is 2-point oscillation in solution)
# one pico = (1000)^4 yocto

rho=0.5
Lambda=2100*umeter**3
Jbeta=0.1e-24*mole/second   # IP3 input
# kappad=1*umolar # <<<<< ???? --> kappa_delta

### Astrocyte parameters
# ---  Calcium fluxes
O_P = 4.4*umolar/second      # Maximal Ca^2+ uptake rate by SERCAs (Op) (0.9)
K_P = 0.05 * umolar          # Ca2+ affinity of SERCAs (Kp)
Omega_C = 6./second           # Maximal rate of Ca^2+ release by IP_3Rs
Omega_L = 0.11/second         # Maximal rate of Ca^2+ leak from the ER (0.1)
# --- IP_3R kinectics
d_1 = 0.13*umolar            # IP_3 binding affinity
d_2 = 1.049*umolar            # Ca^2+ inactivation dissociation constant (1.05)
O_2 = 0.1335/umolar/second      # IP_3R binding rate for Ca^2+ inhibition (0.2)
d_3 = 0.9434*umolar          # IP_3 dissociation constant
d_5 = 0.082*umolar            # Ca^2+ activation dissociation constant (0.08)
# --- IP_3 production
O_delta = 0.02e-24*mole/second  # Maximal rate of IP_3 production by PLCdelta  (0.06)
kappa_delta = 1.0* umolar    # Inhibition constant of PLC_delta by IP_3  (kappad?) (consistent with ode file)
#kappa_delta = 1.5* umolar    # Inhibition constant of PLC_delta by IP_3  (kappad?)
K_delta = 0.1*umolar         # Ca^2+ affinity of PLCdelta
# --- IP_3 degradation
O_5P = 0.05e-24*mole/second       # Maximal rate of IP_3 degradation by IP-5P (O5p)
K_5P = 10.*umolar
K_D = 0.7*umolar             # Ca^2+ affinity of IP3-3K
K_3K = 1.0*umolar            # IP_3 affinity of IP_3-3K
O_3K = 2.0e-24*mole/second     # Maximal rate of IP_3 degradation by IP_3-3K  (4.5)
# --- IP_3 diffusion

#-------------------------------------------------------
# ODEs
#-------------------------------------------------------
# IP3   : cytosolic IP3
# CaCYT : cytosolic calcium
# CaER  : endoplasmic reticulum calcium
# h     : IP3R deinactivation
#-------------------------------------------------------

################################################################################
# Model definition
################################################################################
defaultclock.dt = sim_dt     # Set the integration time

comp_eqs = '''
dIP3/dt=(Jbeta+Jdelta-J3k-J5p)*1e3/(Lambda*(1-rho))  : mmolar
dCaCYT/dt = Jcicr+Jleak-Jserca  : mmolar
dCaER/dt=(-Jcicr-Jleak+Jserca) / (rho / (1-rho) ) : mmolar
dh/dt=(hinf-h)/tauh : 1

# Jdelta=O_delta * ( 1 - ( IP3 /(kappa_delta + IP3) ) * ( CaCYT**2/(CaCYT**2+K_delta**2) )  ) : mmolar/second
# J3k=O_3K * ( ( CaCYT**4/(CaCYT**4+K_D**4)) * (IP3/(K_3K+IP3) ) ) : mmolar/second
# J5p=O_5P * IP3 /(IP3 + K_5P) : mmolar/second
Jdelta=O_delta * ( 1 - ( IP3 /(kappa_delta + IP3) ) * ( CaCYT**2/(CaCYT**2+K_delta**2) )  ) : mole/second
J3k=O_3K * ( ( CaCYT**4/(CaCYT**4+K_D**4)) * (IP3/(K_3K+IP3) ) ) : mole/second
J5p=O_5P * IP3 /(IP3 + K_5P) : mole/second

# Calcium core
minf=IP3/(IP3+d_1) : 1
ninf=CaCYT/(CaCYT+d_5) : 1
Jcicr=Omega_C*(minf*ninf*h)**3*(CaER-CaCYT) : mmolar/second
Jleak=Omega_L*(CaER-CaCYT) : mmolar/second
Jserca=O_P*CaCYT**2/(K_P**2+CaCYT**2) : mmolar/second

# Gating variable
Q2=d_2*(IP3+d_1)/(IP3+d_3) : mmolar
hinf=Q2/(Q2+CaCYT) : 1
tauh=1./O_2/(Q2+CaCYT) : second
'''

#-------------------------------------------------------

N_comp = 1 # Total number of compartments in the network
compartments = NeuronGroup(N_comp, comp_eqs, method='rk4')

# Initial conditions
# How are these set more generally with random numgers? 
compartments.CaCYT = 0.0*umolar
compartments.CaER = 7.5*umolar
compartments.h =  0.95
compartments.IP3 = 0.0*umolar   # Units: umolar or ymolar? 

# # Diffusion between compartments (REVISIT)
# # INCLUDE BEGIN
# comp_to_comp_eqs = '''
# #delta_I = I_post - I_pre        : mmolar
# #J_coupling_post = -F/2 * (1 + tanh((abs(delta_I) - I_Theta)/omega_I)) *
#                   #sign(delta_I) : mmolar/second (summed)
# '''
#
# comp_to_comp = Synapses(compartments, compartments, model=comp_to_comp_eqs)
# # INCLUDE END
# # Couple neighboring compartments: two connections per compartment pair, as
# # the above formulation will only update the I_coupling term of one of the
# # compartments
# # INCLUDE BEGIN
# comp_to_comp.connect('j == (i + 1) % N_pre or '
#                        'j == (i + N_pre - 1) % N_pre')
# # INCLUDE END

################################################################################
# Monitors
################################################################################
astro_mon = StateMonitor(compartments, variables=['CaCYT'], record=True)

################################################################################
# Simulation run
################################################################################
run(duration, report='text')

################################################################################
# Analysis and plotting
################################################################################
#plt.style.use('figures.mplstyle')
# plt.style.use('figures_5.mplstyle')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6.26894, 6.26894 * 0.66),
                       gridspec_kw={'left': 0.1, 'bottom': 0.12})

scaling = 1.2
step = 10

# print(astro_mon.CaCYT)
# print(astro_mon.CaCYT.max(), scaling)
# quit()

# ax.plot(astro_mon.t/second,
#         (astro_mon.CaCYT[0:N_comp//2].T/astro_mon.CaCYT.max() +
#          np.arange(N_comp//2)*scaling), color='black')
# ax.plot(astro_mon.t/second, (astro_mon.CaCYT[N_comp//2+1:].T/astro_mon.CaCYT.max() +
#                              np.arange(N_comp//2+1, N_comp)*scaling),
#         color='black')
# ax.plot(astro_mon.t/second, (astro_mon.CaCYT[N_comp//2].T/astro_mon.CaCYT.max() +
#                              N_comp//2 * scaling),
#         color='C0')
# ax.set(xlim=(0., duration/second), ylim=(0, (N_comp+1.5)*scaling),
#        xticks=np.arange(0., duration/second, 500), xlabel='time (s)',
#        yticks=np.arange(0.5*scaling, (N_comp + 1.5)*scaling, step*scaling),
#        yticklabels=[str(yt) for yt in np.arange(0, N_comp + 1, step)],
#        ylabel='$CaCYT/CaCYT_{max}$ (cell index)')
# pu.adjust_spines(ax, ['left', 'bottom'])
#
# pu.adjust_ylabels([ax], x_offset=-0.08)

#plt.savefig('../text/figures/results/example_5_astro_ring_Figure.eps', dpi=600)  # DELETE
# plt.savefig('evan_single_compartment.pdf', dpi=600)  # DELETE

ax.plot(astro_mon.t/second,astro_mon.CaCYT[0]/umolar)
plt.show()
