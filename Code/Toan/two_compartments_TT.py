#created based on Maurizio's note
import numpy as np
import matplotlib.pylab as plt
from brian2 import *

codegen.target = 'cython'     #use C++ standalone mode

###############################################################################
# Model parameters
###############################################################################
### General parameters
duration = 10**2*second       # Total simulation time
sim_dt = 50*ms               # Integrator/sampling step

### Astrocyte parameters
#-- Geometric constraints
rho_A = 0.18                 # ER-to-cytoplasm volume ratio
rho = rho_A / (1.0 + rho_A)     # ER-to-cytoplasm volume ratio
#-- Calcium fluxes
Omega_C = 6/second              # Maximal rate of calcium release by IP3Rs
Omega_L = 0.11/second           # Maximal rate of calcium leak from the ER
O_P = 2.2*umolar/second         # Maximal rate of calcium uptake by SERCA pumps
K_P = 0.1*umolar                # Calcium affinity of SERCA pumps
#-- Dissociation constants
d_1 = 0.13*umolar               # IP3 binding affinity
d_2 = 1.049*umolar              # Inactivating calcium binding affinity
d_3 = 0.9434*umolar             # IP3 binding affinity (with calcium inactivation)
d_5 = 0.08234*umolar            # Activating calcium binding affinity
O_2 = 0.2/(umolar*second)       # Inactivating calcium binding rate
#-- IP3 production
O_delta = 0.6*umolar/second     # Maximal rate of IP3 production by PLCdelta
kappa_delta = 1.5*umolar        # Inhibition constant of PLCdelta by IP3
K_delta = 0.1*umolar            # Calcium affinity of PLCdelta
#-- IP3 degradation
O_3K = 4.5*umolar/second        # Maximal rate of IP3 degradation by IP3-3K
K_3K = 1.0*umolar               # IP3 affinity of IP3-3K
K_D = 0.7*umolar                # Calcium affinity of IP3-3K
O_5P = 0.05/second              # Maximal rate of IP3 degradation by IP-5P
K_5P = 10*umolar                # ???


################################################################################
# Additional and modified parameters (coppied from Maurizio's code)
################################################################################
# Volume of an average soma
Lambda = 2100*umeter**3

# Multiply value of amplitudes by cytosolic volume (changed rho_A to rho)
O_delta = 0.6  * umolar * Lambda * (1-rho) / second
O_3K    = 4.5  * umolar * Lambda * (1-rho) / second
O_5P    = 0.05 * umolar * Lambda * (1-rho) / second
F       = 0.1 / second
D       = 0.05 / second

print("**** Omega_L= ", Omega_L)
print("**** O_delta= ", O_delta)
print("**** K_delta= ", K_delta)


###############################################################################
# Model definition
###############################################################################

# Astrocytes
astro_eqs = '''
#-- Compartment constraints
rho_A                           : 1  # ratio between ER and compartment volume
rho = rho_A/(1 + rho_A)       : 1  # ratio between ER and cytosolic volume
vol_coef = Lambda * (1 - rho)   : meter**3  # cytosolic volume

#-- Governing equations
dCaCYT/dt = J_chan + J_leak - J_pump                             : mmolar
dCaER/dt = -(J_chan + J_leak - J_pump) / rho_A                   : mmolar
dI/dt = (J_ex + J_delta - J_3K - J_5P + J_coupling) / vol_coef   : mmolar
dh/dt = Omega_h * (h_inf - h)                                    : 1

# where
    J_chan = Omega_C * (m_inf  *  n_inf * h)**3 * (CaER - CaCYT)   : mmolar/second
        m_inf = I/(I + d_1)                                        : 1
        n_inf = CaCYT/(CaCYT + d_5)                                : 1
    J_leak = Omega_L * (CaER - CaCYT)                              : mmolar/second
    J_pump = O_P * CaCYT**2/(CaCYT**2 + K_P**2)                    : mmolar/second
    
    J_ex = -F * delta_I_bias * vol_coef                                        : mole/second
        delta_I_bias = I - stimulus * I_bias                                   : mmolar
            stimulus = int((t % (50*second))<20*second)                        : 1
            I_bias                                                             : mmolar (constant)
    J_delta = O_delta/(1 + I/kappa_delta) * CaCYT**2/(CaCYT**2 + K_delta**2)   : mole/second
    J_3K = O_3K * CaCYT**4/(CaCYT**4 + K_D**4) * I/(I + K_3K)                  : mole/second
    J_5P = O_5P * I/(I + K_5P)                                                 : mole/second
    J_coupling                                                                 : mole/second
    
    h_inf = Q_2/(Q_2 + CaCYT)               : 1                       
    Omega_h = O_2 * (Q_2 + CaCYT)           : Hz
        Q_2 = d_2 * (I + d_1) / (I + d_3)   : mmolar
'''


###############################################################################
# Initialize simulation parameters
###############################################################################
defaultclock.dt = sim_dt     # Set the integration time
N_astro = 2
astrocytes = NeuronGroup(N_astro, astro_eqs, method='rk4')
astrocytes.I_bias = np.asarray([10, 0], dtype=float)*umolar
astrocytes.rho_A = 0.18
astrocytes.CaER = 2*umolar


###############################################################################
# Diffusion between compartments
###############################################################################
astro_to_astro_eqs = '''
delta_I = I_post - I_pre                             : mmolar
J_coupling_post = -D * vol_coef * delta_I            : mole/second (summed)
'''
astro_to_astro = Synapses(astrocytes, astrocytes, model=astro_to_astro_eqs)
astro_to_astro.connect()


###############################################################################
# Prepare data recording
###############################################################################
astro_mon = StateMonitor(astrocytes, variables=['CaCYT', 'CaER', 'J_chan', 'J_leak', 'J_pump',
                                                'I', 'J_ex', 'J_delta', 'J_3K', 'J_5P'], record=True)
    
    
###############################################################################
# Simulation run
###############################################################################
run(duration, report='text')


###############################################################################
# Output and plot data
###############################################################################
fig, axes = plt.subplots(nrows=5, ncols=2, figsize=(2 * 6.26894, 4 * 6.26894 * 0.66),
                       gridspec_kw={'left': 0.1, 'bottom': 0.12})

ax = axes[0, 0]
ax.plot(astro_mon.t/second,astro_mon.CaCYT[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.CaCYT[1]/umolar,'r-')
ax.set_ylabel("CaCYT")
ax = axes[1, 0]
ax.plot(astro_mon.t/second,astro_mon.CaER[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.CaER[1]/umolar,'r-')
ax.set_ylabel("CaER")
ax = axes[2, 0]
ax.plot(astro_mon.t/second,astro_mon.J_chan[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.J_chan[1]/umolar,'r-')
ax.set_ylabel("J_chan")
ax = axes[3, 0]
ax.plot(astro_mon.t/second,astro_mon.J_leak[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.J_leak[1]/umolar,'r-')
ax.set_ylabel("J_leak")
ax = axes[4, 0]
ax.plot(astro_mon.t/second,astro_mon.J_pump[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.J_pump[1]/umolar,'r-')
ax.set_ylabel("J_pump")
ax = axes[0, 1]
ax.plot(astro_mon.t/second,astro_mon.I[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.I[1]/umolar,'r-')
ax.set_ylabel("I")
ax = axes[1, 1]
ax.plot(astro_mon.t/second,astro_mon.J_ex[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.J_ex[1]/umolar,'r-')
ax.set_ylabel("J_ex")
ax = axes[2, 1]
ax.plot(astro_mon.t/second,astro_mon.J_delta[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.J_delta[1]/umolar,'r-')
ax.set_ylabel("J_delta")
ax = axes[3, 1]
ax.plot(astro_mon.t/second,astro_mon.J_3K[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.J_3K[1]/umolar,'r-')
ax.set_ylabel("J_3K")
ax = axes[4, 1]
ax.plot(astro_mon.t/second,astro_mon.J_5P[0]/umolar,'b-')
ax.plot(astro_mon.t/second,astro_mon.J_5P[1]/umolar,'r-')
ax.set_ylabel("J_5P")

plt.tight_layout()
plt.savefig("plot.pdf", bbox_inches='tight')