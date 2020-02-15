# 2020-02-04,1.15pm: Started synapse_groups branch. 
# Objective: refine the two-synapse group architecture
# 2020-02-07: add IP3 and Ce equations, even though the code is not certified to be correct. 

#  SIMPLIFY THE CODE TO ONLY KEEP THE Calcium and look at diffusion term only. 

#from brian2 import *
import utils as u
#codegen.target = 'cython'

# Author: Gordon Erlebacher
# Objective: create a pure python code with two compartments as a way to test the Brian2 
# code, which is not stable. 
#

### General parameters
# if dt=1*ms, I get Nan
# VVT is very non-smooth as a function of time
# Following Maurizio's derivation, implement the solution to the cubic equation in Cytosolic Calcium concentration.
# Otherwise very similar to voltage_cubic.py

# Initial conditions should be self-consistent in some way. 
# What is V0, C0, etc at t=0?

class Parameters:
    def __init__(self):
        self.duration = 0.7 #*ms          # Total simulation time
        self.sim_dt = 0.1 #*ms            # Integrator/sampling step
        
        # Model definition
        #self.defaultclock.dt = self.sim_dt     # Set the integration time
        
        # Extracellular calcium (for Goldman formula). 
        # What is the location? Where was Goldman's formula applied? Across the extracelluar membrane (radius R), 
        # or between the cytosol and the ER? 
        self.Ce  = 1.e-4  # * umole / cm**3
        print("Ce= ", self.Ce)
        self.Ce  = 1800  #* uM   # too large? 
        print("Ce= ", self.Ce)
        # From Oschmann Master thesis, page 53. 
        self.Crest = 0.073 # * uM  # M is molar, u is micro
        print("Crest= ", self.Crest)
        
        self.P_Ca = 4.46e-13 # * cm / second   # From "Calcium Permeability in Large Unilamellar Vesicles Prepared from Bovine Lens Cortical Lipids"
        self.V_T = 26 # * mvolt
        self.Vrest  = -80 # * mvolt    # in  your derivation, you set Vrest = 0
        self.F = 96485.3329 # * amp * second / mole  #  is Coulomb
        self.k_B = 1.38064852e-23 # * meter**2 * kilogram / second**2 / kelvin
        self.N_A = 6.02214076e23 # / mole
        self.e = 1.602176634e-19 # * amp * second
        self.Cm = 1 # * uF / meter**2 #: 1  #Capacitance per meter**2
        self.D_C  = 5.3e-6 # * cm**2 / second # from "Free diffusion coefficient of ionic calcium in cytoplasm", Donahue et al (1987)
        self.D_CE = 5.3e-6 # * cm**2 / second # (a guess, ER)
        self.D_I  = 5.3e-6 #* cm**2 / second # (a guess, IP3)
        self.Rgas = 8.31 # * joule / kelvin / mole
        
        
        # Constants found in currents
        # IP3 production
        self.O_delta = 0.6 # *umolar/second  # Maximal rate of IP_3 production by PLCdelta
        self.k_delta = 1.5 #* umolar        # Inhibition constant of PLC_delta by IP_3
        self.K_delta = 0.1 #*umolar         # Ca^2+ affinity of PLCdelta   <<<<<<<<<<<<<<<<<
        
        # IP3 degradation
        self.K_D  = 0.7 # *umolar            # Ca^2+ affinity of IP3-3K
        self.K_3  = 1.0 # *umolar            # IP_3 affinity of IP_3-3K
        self.O_3K = 4.5 # *umolar/second     # Maximal rate of IP_3 degradation by IP_3-3K
        
        # IP5 degradation
        self.o_5P = 0.05 #/second           # Maximal rate of IP_3 degradation by IP-5P
        self.K_5P = 10 #*umolar
        
        # IP3 delta Production
        self.o_delta = 0.6 #*umolar/second  # Maximal rate of IP_3 production by PLCdelta
        self.k_delta = 1.5 #* umolar    # Inhibition constant of PLC_delta by IP_3
        self.K_delta = 0.1 #*umolar         # Ca^2+ affinity of PLCdelta   <<<<<<<<<<<<<<<<<
        
        # Volume of an average soma
        self.Lambda = 2100 #*umeter**3
        # Not sure about this
        self.rho_A = 0.18                 # ER-to-cytoplasm volume ratio ?
        self.rho = self.rho_A / (1.+self.rho_A)     # ER-to-cytoplasm volume ratio ?
        
        # Multiply value of amplitudes by cytosolic volume (changed rho_A to rho)
        self.o_delta = 0.6  # * umolar * Lambda * (1-rho) / second
        self.o_3K    = 4.5  # * umolar * Lambda * (1-rho) / second
        self.o_5P    = 0.05 # * umolar * Lambda * (1-rho) / second
        
        self.Omega_L = 0.1 #/second         # Maximal rate of Ca^2+ leak from the ER
        self.Omega_2 = 0.2 #/umolar/second      # IP_3R binding rate for Ca^2+ inhibition
        self.Omega_u = 0.2 #/umolar/second      # uptake/dissassociation constant
        
        # --- IP_3R kinectics
        self.d_1 = 0.13 #*umolar            # IP_3 binding affinity
        self.d_2 = 1.05 #*umolar            # Ca^2+ inactivation dissociation constant
        self.O_2 = 0.2 #/umolar/second      # IP_3R binding rate for Ca^2+ inhibition
        self.d_3 = 0.9434 #*umolar          # IP_3 dissociation constant
        self.d_5 = 0.08 #*umolar            # Ca^2+ activation dissociation constant
    
        self.p_open = 1
        self.P_r     = 1 #*umolar/second
        self.P_CE    = 1 #*umolar/second
    
        self.eta_p = 1
        self.d_ER  = 1  #* umeter
    
        # Another complex calculation that depends on solution to a cubic equation
        # eqs. (73)-(76) in De Pitta's notes. MUST IMPLEMNENT. 
        # Later, evaluate the terms and figure out what is negligeable and what is not. 
        self.dv_ER = 0  #* mvolt
    
        self.K_P = 0.05  #* umolar          # Ca2+ affinity of SERCAs
    
        #o_3K
        #p_open
        #Pr
        #dv_ER
        #d_ER
        #N_A
        #Omega_u
        #eta_p

        ################################################################################
        # Additional and modified parameters
        ################################################################################
        # Volume of an average soma
        self.Lambda = 2100   #*umeter**3
        
        # Multiply value of amplitudes by cytosolic volume (changed rho_A to rho)
        self.O_delta = 0.6  #  * umolar * Lambda * (1-rho) / second
        self.O_3K    = 4.5  #  * umolar * Lambda * (1-rho) / second
        self.O_5P    = 0.05 #  * umolar * Lambda * (1-rho) / second
        self.F       = 0.1  #  / second
        self.D       = 0.05 #  / second  # setting D to zero has no effect. Why?
        
        print("**** Omega_L= ", self.Omega_L)
        print("**** O_delta= ", self.O_delta) # 0.6
        print("**** K_delta= ", self.K_delta) # 0.6
#----------------------------------------------------------------------
class RHS(): 
    # Right-hand side of differential equations
    # 

    def __init__(self, g):
        self.g = g

    def rhs(self, C, Ce, h, I):
        # C: calcium
        # Ce: calcium in ER
        # h : fraction of open channels
        # h : IP3 concentration

        g = self.g
        J_C = diff_C = 0
        J_Ce = diff_Ce = 0
        J_I = diff_I = 0

        rhs_C  = J_C + diff_C
        rhs_Ce = J_Ce + diff_Ce
        rhs_I  = J_I + diff_I

        #dI/dt = (Jbeta + Jdelta - J3K - J5P) / (Lambda*(1-rho)) + 1*coupling_CE   : mole/meter**3
        #dC/dt = 1*coupling_C + 1.*coupling_electro + 0.*electro_diffusion + Jr + J1 - Jp  : mole / meter**3
        #dCE/dt = 1*coupling_CE + Jp/(rho*Lambda) - (Jr + J1)  : mole/meter**3

        #IP3 dynamics
        Jbeta  = 0  #*mole/meter**3/second  : mole/meter**3/second
        Jdelta = g.o_delta * (g.k_delta/(I+g.k_delta)) * (C**2/(C**2+g.K_delta**2)) #: mole/second  # not sure about units, o_delta, k_delta, K_delta
        J5P    = g.o_5P * (I/(I+g.K_5P)) #: mole/second # o_5P, K_5P
        J3K    = g.o_3K * (C**4/(C**4+g.K_D**4)) * (I/(I+g.K_3)) #: mole/second # o_3K, K_D, K_3

        coupling_Ce = 0 # CHECK
        rhs_I = (Jbeta + Jdelta - J3K - J5P) / (g.Lambda*(1-g.rho)) + 1*coupling_Ce   #: mole/meter**3

        OmegaH = (g.Omega_2*(I+g.d_1) + g.O_2*(I+g.d_3)*C) / (I + g.d_3) #: Hz # Omega_2, O_2, d_1, d_3
        hinf   =  g.d_2 * (I + g.d_1) / (g.d_2*(I + g.d_1) + (I + g.d_3)*C) #: 1 # d_2, d_1, d_3
        rhs_h = OmegaH * (hinf - h) #: 1
#----------------------------------------------------------------------
"""
# has large impact on diffusion. Creates instability
coupling_electro   : mole / second / meter**3
#coupling_electro = 1  : 1 

# Tot_C is the total calcium (computed in synapses2)
Tot_C : mole/meter**3

V = Vrest + (F*s/Cm) * (Tot_C-Crest) : volt
VVT = -2.*V/V_T : 1  # value of 0, which leads to a zero denominator

# Has minimal impact on Calcium
electro_diffusion = -P_Ca * V / (R * V_T) * (Ce*exp(-2*V/V_T) - Tot_C) / (exp(-2*V/V_T) - 1) : mole / second / meter**3

A1 : meter**4 / mole 
B1 : meter 
C1 : mole / meter**2 

CC0 : mole / meter**3
V0 : volt
nb_connections : 1  # number of synaptic connections

# Calcium diffusion
# Each term will handle half a compartment 

dC/dt = 1*coupling_C + 1.*coupling_electro + 0.*electro_diffusion + Jr + J1 - Jp  : mole / meter**3

####  ENDOPLASMIC RETICULUM
Tot_CE                 : mole/meter**3
coupling_CE            : mole/second/meter**3
dCE/dt = 1*coupling_CE + Jp/(rho*Lambda) - (Jr + J1)  : mole/meter**3

####  IP3
Tot_I                  : mole/meter**3
coupling_I             : mole/second/meter**3
dI/dt = (Jbeta + Jdelta - J3K - J5P) / (Lambda*(1-rho)) + 1*coupling_CE   : mole/meter**3

### Open Channels
dh/dt = OmegaH * (hinf - h) : 1


################3
###  CURRENTS

Jr     = (2*r/dR2) * P_r * p_open * (CE-C) : mole/meter**3/second  # p_open, Pr
J1     = (4*P_CE/r*VVT) * dv_ER * (C*exp(-2.*dv_ER/VVT) - CE) / (exp(-2.*dv_ER/VVT)-1.) : mole/meter**3/second  # dv_ER
Jp     = (2*r*d_ER)/(N_A*dR2) * Omega_u * eta_p * C**2 / (C**2 + K_P**2) : mole/meter**3/second # d_ER, N_A, Omega_u, eta_p, K_p
minf   =  I / (I + d_1)  : 1 # d_1
ninf   = C / (C + d_5) : 1 # d_5
hinf   =  d_2 * (I + d_1) / (d_2*(I +d_1) + (I+d_3)*C) : 1 # d_2, d_1, d_3
OmegaH = (Omega_2*(I+d_1) + O_2*(I+d_3)*C) / (I + d_3) : Hz # Omega_2, O_2, d_1, d_3



#IP3 dynamics
Jbeta  = 0*mole/meter**3/second  : mole/meter**3/second
Jdelta = o_delta * (k_delta/(I+k_delta)) * (C**2/(C**2+K_delta**2)) : mole/second  # not sure about units, o_delta, k_delta, K_delta
J5P    = o_5P * (I/(I+K_5P)) : mole/second # o_5P, K_5P
J3K    = o_3K * (C**4/(C**4+K_D**4)) * (I/(I+K_3)) : mole/second # o_3K, K_D, K_3

'''

N_astro = 3 # Total number of astrocytes1 in the network (each compartment broken into two
N_astro = 2*N_astro # Total number of astrocytes1 in the network (each compartment broken into two
astrocytes1 = NeuronGroup(N_astro, astro_eqs, method='euler', name='ng1', order=1)

# Initital Conditions
astrocytes1.L = 8 *  umeter  # length of a compartment
astrocytes1.R = 3 * umeter
astrocytes1.r = 2 * umeter
# Calcium in 1/2 compartment should be the value for the full compartment since the 1/2 compartment
# Each 1/2 compartment only considers the fork at one end. 
# is an artifact of the numerical scheme.
astrocytes1.C = 1.1e-4  * mole / meter**3
astrocytes1.I = 1.1e-4  * mole / meter**3
astrocytes1.CE = 1.1e-4  * mole / meter**3
#astrocytes1.C = [1.1e-4, 1.5e-4, 1.6e-4]*2 * mole / meter**3

astro_mon = StateMonitor(astrocytes1, variables=['Tot_C', 'Tot_CE', 'Tot_I', 'coupling_C', 'coupling_electro', 'electro_diffusion', 'nb_connections', 'A1', 'B1', 'C1', 'CC0', 'V0'], record=True)
b_mon = StateMonitor(astrocytes1, variables=['CC0', 'dR2', 's', 'L'], record=True)

'''
One issue is Boundary Conditions. I have not done anything to implement those. 
Also, each synapse has its own value of V0. 
'''

#----------------------------------------------------------------------
# Diffusion between astrocytes1
synapse1_eqs = '''
# Assuming the above is correct, let us figure out units. 
# using C_pre:  C * F * s / Cm = V  ==> F / Cm = V / s / C = (V / s) * m^3 / mole

####  CYTOPLASM

# A1 stable: value does not change from iteration to iteration
A1_pre = (0.5 * F * s_post)/(Cm * V_T) * (dR2_post / L_post) : meter**4 / mole (summed)
B1_pre = (1. - (s_post * F * Tot_C_post)/(2. * Cm * V_T)) : meter (summed)
C1_pre = dR2_post * Tot_C_post / L_post : mole / meter**2 (summed)

# Only computed to trick Brian into computing C0_pre by summing the values for all synapses (which will 
#   all bequal  since it depends on A1, B1, C1, for which the summation has already been executed. 
"""
#----------------------------------------------------------------------
g = Parameters()
rhs = RHS(g)
rhs.rhs(1.,1.,1.,1.)


#----------------------------------------------------------------------
astro_eqs = '''
s = R + r          : meter
r                  : meter
dR2 = R**2 - r**2  : meter**2
DR2                : meter**2
L                  : meter
R                  : meter  # Outer radius
coupling_C         : mole/second/meter**3

####  CALCIUM

# has large impact on diffusion. Creates instability
coupling_electro   : mole / second / meter**3
#coupling_electro = 1  : 1 

# Tot_C is the total calcium (computed in synapses2)
Tot_C : mole/meter**3

V = Vrest + (F*s/Cm) * (Tot_C-Crest) : volt
VVT = -2.*V/V_T : 1  # value of 0, which leads to a zero denominator

# Has minimal impact on Calcium
electro_diffusion = -P_Ca * V / (R * V_T) * (Ce*exp(-2*V/V_T) - Tot_C) / (exp(-2*V/V_T) - 1) : mole / second / meter**3

A1 : meter**4 / mole 
B1 : meter 
C1 : mole / meter**2 

CC0 : mole / meter**3
V0 : volt
nb_connections : 1  # number of synaptic connections

# Calcium diffusion
# Each term will handle half a compartment 

dC/dt = 1*coupling_C + 1.*coupling_electro + 0.*electro_diffusion + Jr + J1 - Jp  : mole / meter**3

####  ENDOPLASMIC RETICULUM
Tot_CE                 : mole/meter**3
coupling_CE            : mole/second/meter**3
dCE/dt = 1*coupling_CE + Jp/(rho*Lambda) - (Jr + J1)  : mole/meter**3

####  IP3
Tot_I                  : mole/meter**3
coupling_I             : mole/second/meter**3
dI/dt = (Jbeta + Jdelta - J3K - J5P) / (Lambda*(1-rho)) + 1*coupling_CE   : mole/meter**3

### Open Channels
dh/dt = OmegaH * (hinf - h) : 1


################3
###  CURRENTS

Jr     = (2*r/dR2) * P_r * p_open * (CE-C) : mole/meter**3/second  # p_open, Pr
J1     = (4*P_CE/r*VVT) * dv_ER * (C*exp(-2.*dv_ER/VVT) - CE) / (exp(-2.*dv_ER/VVT)-1.) : mole/meter**3/second  # dv_ER
Jp     = (2*r*d_ER)/(N_A*dR2) * Omega_u * eta_p * C**2 / (C**2 + K_P**2) : mole/meter**3/second # d_ER, N_A, Omega_u, eta_p, K_p
minf   =  I / (I + d_1)  : 1 # d_1
ninf   = C / (C + d_5) : 1 # d_5
hinf   =  d_2 * (I + d_1) / (d_2*(I +d_1) + (I+d_3)*C) : 1 # d_2, d_1, d_3
OmegaH = (Omega_2*(I+d_1) + O_2*(I+d_3)*C) / (I + d_3) : Hz # Omega_2, O_2, d_1, d_3



#IP3 dynamics
Jbeta  = 0*mole/meter**3/second  : mole/meter**3/second
Jdelta = o_delta * (k_delta/(I+k_delta)) * (C**2/(C**2+K_delta**2)) : mole/second  # not sure about units, o_delta, k_delta, K_delta
J5P    = o_5P * (I/(I+K_5P)) : mole/second # o_5P, K_5P
J3K    = o_3K * (C**4/(C**4+K_D**4)) * (I/(I+K_3)) : mole/second # o_3K, K_D, K_3

'''

N_astro = 3 # Total number of astrocytes1 in the network (each compartment broken into two
N_astro = 2*N_astro # Total number of astrocytes1 in the network (each compartment broken into two
astrocytes1 = NeuronGroup(N_astro, astro_eqs, method='euler', name='ng1', order=1)

# Initital Conditions
astrocytes1.L = 8 *  umeter  # length of a compartment
astrocytes1.R = 3 * umeter
astrocytes1.r = 2 * umeter
# Calcium in 1/2 compartment should be the value for the full compartment since the 1/2 compartment
# Each 1/2 compartment only considers the fork at one end. 
# is an artifact of the numerical scheme.
astrocytes1.C = 1.1e-4  * mole / meter**3
astrocytes1.I = 1.1e-4  * mole / meter**3
astrocytes1.CE = 1.1e-4  * mole / meter**3
#astrocytes1.C = [1.1e-4, 1.5e-4, 1.6e-4]*2 * mole / meter**3

astro_mon = StateMonitor(astrocytes1, variables=['Tot_C', 'Tot_CE', 'Tot_I', 'coupling_C', 'coupling_electro', 'electro_diffusion', 'nb_connections', 'A1', 'B1', 'C1', 'CC0', 'V0'], record=True)
b_mon = StateMonitor(astrocytes1, variables=['CC0', 'dR2', 's', 'L'], record=True)

'''
One issue is Boundary Conditions. I have not done anything to implement those. 
Also, each synapse has its own value of V0. 
'''

#----------------------------------------------------------------------
# Diffusion between astrocytes1
synapse1_eqs = '''
# Assuming the above is correct, let us figure out units. 
# using C_pre:  C * F * s / Cm = V  ==> F / Cm = V / s / C = (V / s) * m^3 / mole

####  CYTOPLASM

# A1 stable: value does not change from iteration to iteration
A1_pre = (0.5 * F * s_post)/(Cm * V_T) * (dR2_post / L_post) : meter**4 / mole (summed)
B1_pre = (1. - (s_post * F * Tot_C_post)/(2. * Cm * V_T)) : meter (summed)
C1_pre = dR2_post * Tot_C_post / L_post : mole / meter**2 (summed)

# Only computed to trick Brian into computing C0_pre by summing the values for all synapses (which will 
#   all bequal  since it depends on A1, B1, C1, for which the summation has already been executed. 
# There must be an easier method without calculating C0 and V0 in the NeuronGroup
nb_connections_pre = 1 : 1 (summed)

# Notice I am dividing by the number of connections (3 if a single fork)
CC0_pre = ((-B1_pre + sqrt(B1_pre**2 + 4*C1_pre*A1_pre)) / (2.*A1_pre)) / nb_connections_pre : mole / meter**3 (summed)
V0_pre = (Vrest + (CC0_pre - Crest) * (F * s) / Cm) / nb_connections_pre   : volt  (summed)

# Investigate why this should be a minus sign
coupling_C_post = (4*D_C / L_post**2) * (Tot_C_pre - Tot_C_post) : mole/second/meter**3 (summed)

# MAKE SURE CC0 and V0 are computed before updating the coupling parameters. 
coupling_electro_pre = (4*D_C/L_post**2/V_T) * (CC0_post + Tot_C_post) * (V0_post - V_post) : mole/second/meter**3 (summed) 

####  ENDOPLASMIC RETICULUM
coupling_CE_post = (4*D_CE / L_post**2) * (Tot_CE_pre - Tot_CE_post) : mole/second/meter**3 (summed)

####  IP3
coupling_I_post = (4*D_I / L_post**2) * (Tot_I_pre - Tot_I_post) : mole/second/meter**3 (summed)
'''
#----------------------------------------------------------------------

# TEMPORARY
synapse2_eqs = '''
    # C is updated in each compartment of the astrocyte and combined here
    Tot_C_syn   = C_pre  + C_post    : mole / meter**3
    Tot_C_pre   = Tot_C_syn       : mole / meter**3 (summed)

    # update Calcium in compartments
    # This works because synapses are bidirectional i --> j and j --> i

####  ENDOPLASMIC RETICULUM
    Tot_CE_pre  = Tot_CE_syn      : mole / meter**3 (summed)
    Tot_CE_syn  = CE_pre + CE_post   : mole / meter**3

####  IP3
    Tot_I_syn = I_pre + I_post : mole / meter**3
    Tot_I_pre = Tot_I_syn      : mole / meter**3 (summed)
'''

synapses1 = Synapses(astrocytes1, astrocytes1, model=synapse1_eqs, method='euler', order=0, name='sg1')
synapses2 = Synapses(astrocytes1, astrocytes1, model=synapse2_eqs, method='euler', order=5, name='sg2')

# Connections count from 0
#  ---*---.---*----
#         |___*____

# Single fork (N_astro = 6)
# Six Compartments: 0,1,...,4,5
for i in range(N_astro):
    synapses1.connect(i=i,j=i)
pairs = [(3,1),(3,2),(1,2)]
for pair in pairs:
    synapses1.connect(i=pair[0], j=pair[1])
    synapses1.connect(i=pair[1], j=pair[0])

pairs = [(0,1),(2,4),(3,5)]
for pair in pairs:
    synapses2.connect(i=pair[0], j=pair[1])
    synapses2.connect(i=pair[1], j=pair[0])

s = synapses1.summed_updaters
s['nb_connections_pre'].order = 0
s['A1_pre'].order = 2
s['B1_pre'].order = 2
s['C1_pre'].order = 2
s['CC0_pre'].order = 4
s['V0_pre'].order = 5
s['coupling_C_post'].order = 6
s['coupling_electro_pre'].order = 7

astrocytes1.state_updater.order = 10

s = synapses2.summed_updaters
s['Tot_C_pre'].order = 13


for k,v in s.items():
    print("updater[%s]= " % k, v)


# Must call StateMonitor AFTER synaptic connections are established
syn1_mon = StateMonitor(synapses1, variables=[], record=True, name='syn1_mon')
syn2_mon = StateMonitor(synapses2, variables=[], record=True, name='syn2_mon')

groups = [synapses1, synapses2]  #  must have i,j connections
u.connectionMatrices(groups)


#-------------------------------------------
# Initial Conditions (call after connect())
# Temporary. Should be computed in the synapse
#synapses1.P                  = 1.1e-2 

groups = [astrocytes1, synapses1, synapses2]
u.lowLevelInfo(groups)

# Run Simulation
print(scheduling_summary())
run(duration)
#run(duration, report='text')

u.printData(astro_mon, b_mon, syn1_mon, syn2_mon)
quit()
#u.plots(astro_mon)
#----------------------------------------------------------------------

