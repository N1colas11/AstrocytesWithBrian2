# 2020-02-04,1.15pm: Started synapse_groups branch.
# Objective: refine the two-synapse group architecture
# 2020-02-07: add IP3 and Ce equations, even though the code is not certified to be correct.

#  SIMPLIFY THE CODE TO ONLY KEEP THE Calcium and look at diffusion term only.

#from brian2 import *
import utils as u
import numpy as np
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
        self.Ce  = 1.e-4  # * umole / cm**3, ==> MKS:   10^(-6)/ 10^(-6) = mole / meter**3
        print("Ce= ", self.Ce)
        #self.Ce  = 1800  #* uM   # too large?  = 10^-6 M  (what is M versus mole?) <<<<<<<
        self.Ce  = 1.800  #  mole/meter**3
        print("Ce= ", self.Ce)

        #b1 = 1*uM            # = 1 uM
        #c1 = 1*mole / meter**3
        #print("c1/b1= ", c1/b1, "  (1*mole/meter**3) / uM") # 1*uM = 10^{-3} * mole/meter**3

        # From Oschmann Master thesis, page 53.
        #self.Crest = 0.073 # * uM  # M is molar, u is micro = 10^(-6) M  (M: Molar)
        self.Crest = 0.073*1.e-3 # * uM  # M is molar, u is micro = 10^(-6) M  (M: Molar)

        #self.P_Ca = 4.46e-13 # * cm / second = 10^-2 m*Hz # From "Calcium Permeability in Large Unilamellar Vesicles Prepared from Bovine Lens Cortical Lipids"
        self.P_Ca = 4.46e-15 # * m / second

        #self.V_T = 26 # * mvolt = 1.e-3 V
        self.V_T = 0.026 # * volt

        #self.Vrest  = -80 # * mvolt    # in  your derivation, you set Vrest = 0
        self.Vrest  = -80e-3 # * volt    # in  your derivation, you set Vrest = 0

        self.F = 96485.3329 # * amp * second / mole  #  is Coulomb
        self.k_B = 1.38064852e-23 # * meter**2 * kilogram / second**2 / kelvin
        self.N_A = 6.02214076e23 # / mole
        self.e = 1.602176634e-19 # * amp * second

        #self.Cm = 1 # * uF / meter**2 #: 1  #Capacitance per meter**2 = 10^-6 F / meter**2
        self.Cm = 1e-6 # * uF / meter**2 #: 1  #Capacitance per meter**2 = 10^-6 F / meter**2

        #self.D_C  = 5.3e-6 # * cm**2 / second # from "Free diffusion coefficient of ionic calcium in cytoplasm", Donahue et al (1987)
        self.D_C  = 5.3e-10 # * m**2 / second # from "Free diffusion coefficient of ionic calcium in cytoplasm", Donahue et al (1987)

        #self.D_CE = 5.3e-6 # * cm**2 / second # (a guess, ER) = 10-4 meter**2 / second
        self.D_CE = 5.3e-20 # * cm**2 / second # (a guess, ER) = 10-4 meter**2 / second

        #self.D_I  = 5.3e-6 #* cm**2 / second # (a guess, IP3)
        self.D_I  = 5.3e-10 #* cm**2 / second # (a guess, IP3)

        self.Rgas = 8.31 # * joule / kelvin / mole


        # Constants found in currents
        # IP3 production
        #self.O_delta = 0.6 # *umolar/second  # Maximal rate of IP_3 production by PLCdelta
        self.O_delta = 0.6e-6 # molar/second  # Maximal rate of IP_3 production by PLCdelta

        #self.k_delta = 1.5 #* umolar        # Inhibition constant of PLC_delta by IP_3
        self.k_delta = 1.5e-3 #* mole/meter**3        # Inhibition constant of PLC_delta by IP_3

        #self.K_delta = 0.1 #*umolar         # Ca^2+ affinity of PLCdelta   <<<<<<<<<<<<<<<<<
        self.K_delta = 0.1e-3 #*mole/meter**3         # Ca^2+ affinity of PLCdelta   <<<<<<<<<<<<<<<<<

        # IP3 degradation
        #self.K_D  = 0.7 # *umolar            # Ca^2+ affinity of IP3-3K
        self.K_D  = 0.7e-3 # *mole/meter**3            # Ca^2+ affinity of IP3-3K

        #self.K_3  = 1.0 # *uM            # IP_3 affinity of IP_3-3K
        self.K_3  = 1.0e-3 # *mole/meter**3            # IP_3 affinity of IP_3-3K

        #self.O_3K = 4.5 # *umolar/second     # Maximal rate of IP_3 degradation by IP_3-3K
        self.O_3K = 4.5e-3 # mole/meter^3/second     # Maximal rate of IP_3 degradation by IP_3-3K

        # IP5 degradation
        self.o_5P = 0.05 #/second           # Maximal rate of IP_3 degradation by IP-5P
        #self.K_5P = 10 #*umolar
        self.K_5P = 10e-3 #*mole/meter**3

        # IP3 delta Production
        #self.o_delta = 0.6 #*umolar/second  # Maximal rate of IP_3 production by PLCdelta
        self.o_delta = 0.6e-3 #*mole/meter**3/second  # Maximal rate of IP_3 production by PLCdelta
        #self.k_delta = 1.5 #* umolar    # Inhibition constant of PLC_delta by IP_3
        self.k_delta = 1.5e-3 #* mole/meter**3/second    # Inhibition constant of PLC_delta by IP_3
        #self.K_delta = 0.1 #*umolar         # Ca^2+ affinity of PLCdelta   <<<<<<<<<<<<<<<<<
        self.K_delta = 0.1e-3 #*mole/meter**3         # Ca^2+ affinity of PLCdelta   <<<<<<<<<<<<<<<<<
        #print("1uM/1umolar= ", 1*uM/(1*umolar)); # answer is 1

        # Volume of an average soma
        #self.Lambda = 2100 #*umeter**3
        self.Lambda = 2100e-18 #meter**3
        # Not sure about this
        self.rho_A = 0.18                 # ER-to-cytoplasm volume ratio ?
        self.rho = self.rho_A / (1.+self.rho_A)     # ER-to-cytoplasm volume ratio ?

        # Multiply value of amplitudes by cytosolic volume (changed rho_A to rho)
        # It does not look as if I did the multiplication
        #self.o_delta = 0.6  # * umolar * Lambda * (1-rho) / second
        self.o_delta = 0.6e-3  # * mole/meter**3 * Lambda * (1-rho) / second
        #self.o_3K    = 4.5  # * umolar * Lambda * (1-rho) / second
        self.o_3K    = 4.5e-3  # * mole/meter**3 * Lambda * (1-rho) / second
        self.o_5P    = 0.05e-3 # * mole/meter**3 * Lambda * (1-rho) / second

        self.Omega_L = 0.1   #/second         # Maximal rate of Ca^2+ leak from the ER
        self.Omega_2 = 0.2   #/second      # IP_3R binding rate for Ca^2+ inhibition
        self.Omega_u = 0.2-3 #/second      # uptake/dissassociation constant

        # --- IP_3R kinectics
        #self.d_1 = 0.13 #*umolar            # IP_3 binding affinity
        self.d_1 = 0.13e-3 #*mole/meter**3            # IP_3 binding affinity
        #self.d_2 = 1.05 #*umolar            # Ca^2+ inactivation dissociation constant
        self.d_2 = 1.05e-3 #*mole/meter**3            # Ca^2+ inactivation dissociation constant
        #self.O_2 = 0.2 #/umolar/second      # IP_3R binding rate for Ca^2+ inhibition
        self.O_2 = 0.2e3 #/(mole/meter**3)/second      # IP_3R binding rate for Ca^2+ inhibition
        #self.d_3 = 0.9434 #*umolar          # IP_3 dissociation constant
        self.d_3 = 0.9434e-3 #*mole/meter**3          # IP_3 dissociation constant
        #self.d_5 = 0.08 #*umolar            # Ca^2+ activation dissociation constant
        self.d_5 = 0.08 #*mole/meter**3            # Ca^2+ activation dissociation constant

        self.p_open = 1
        #self.P_r     = 1 #*umolar/second
        self.P_r     = 1e-3 #*mole/meter**3/second
        #self.P_CE    = 1 #*umolar/second
        self.P_CE    = 1e-3 #*mole/meter**3/second

        self.eta_p = 1
        #self.d_ER  = 1  #* umeter
        self.d_ER  = 1e-6  #* meter

        # Another complex calculation that depends on solution to a cubic equation
        # eqs. (73)-(76) in De Pitta's notes. MUST IMPLEMNENT.
        # Later, evaluate the terms and figure out what is negligeable and what is not.
        #self.dv_ER = 10  #* mvolt  (CANNOT BE zero, else a NaN will result)
        self.dv_ER = 10e-3  #* mvolt  (CANNOT BE zero, else a NaN will result)

        #self.K_P = 0.05  #* umolar          # Ca2+ affinity of SERCAs
        self.K_P = 0.05e-3  #* mole/meter**3/          # Ca2+ affinity of SERCAs

        # variables normally vary per compartment
        #astrocytes1.L = 8 *  umeter  # length of a compartment
        #astrocytes1.R = 3 * umeter
        #astrocytes1.r = 2 * umeter
        #self.L = 8 # *umeter
        self.L = 8e-6 # *meter
        #self.R = 3 # *umeter
        self.R = 3e-6 # *umeter
        #self.r = 2 # *umeter
        self.r = 2e-6 # *umeter

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
        #self.Lambda = 2100   #*umeter**3
        self.Lambda = 2100e-18   #*meter**3

        # Multiply value of amplitudes by cytosolic volume (changed rho_A to rho)
        #self.O_delta = 0.6  #  * umolar * Lambda * (1-rho) / second
        self.O_delta = 0.6e-21  #  * mole/meter**3 * Lambda * (1-rho) / second
        #self.O_3K   = 4.5  #   * umolar * Lambda * (1-rho) / second
        self.O_3K    = 4.5e-21  #   * mole/meter**3 * Lambda * (1-rho) / second
        #self.O_5P    = 0.05 #  * umolar * Lambda * (1-rho) / second
        self.O_5P    = 0.05e-21 #  * mole/meter**3 * Lambda * (1-rho) / second
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
        g = self.g

        # C: calcium
        # Ce: calcium in ER
        # h : fraction of open channels
        # h : IP3 concentration

        # rhs_C, etc. Are arrays, created via implicit vectorization involving numpy arrays

        dR2 = g.R**2 - g.r**2  #: meter**2
        s = g.R + g.r          #: meter

        V = g.Vrest + (g.F*s/g.Cm) * (C-g.Crest) #: volt
        VVT = -2.*V/g.V_T #: 1  # value of 0, which leads to a zero denominator

        J_C = diff_C = 0
        J_Ce = diff_Ce = 0
        J_I = diff_I = 0

        rhs_C  = J_C + diff_C
        rhs_Ce = J_Ce + diff_Ce
        rhs_I  = J_I + diff_I

        #dI/dt = (Jbeta + Jdelta - J3K - J5P) / (Lambda*(1-rho)) + 1*coupling_CE   : mole/meter**3
        #dC/dt = 1*coupling_C + 1.*coupling_electro + 0.*electro_diffusion + Jr + J1 - Jp  : mole / meter**3
        #dCE/dt = 1*coupling_CE + Jp/(rho*Lambda) - (Jr + J1)  : mole/meter**3

        # Cytosol Dynamics
        Tot_C = C
        Tot_CE = Ce               #: mole/meter**3  (sum of Ce in two half compartments in Brian2 code)

        Jr     = (2*g.r/dR2) * g.P_r * g.p_open * (Ce-C) #: mole/meter**3/second  # p_open, Pr
        #J1     = (4*g.P_CE/g.r*VVT) * g.dv_ER * (C*np.exp(-2.*g.dv_ER/VVT) - Ce) / (np.exp(-2.*g.dv_ER/VVT)-1.) #: mole/meter**3/second  # dv_ER
        print(g.dv_ER, ", VVT= ", VVT)
        xxx = -2.*g.dv_ER/VVT
        J1     = (4*g.P_CE/g.r*VVT) * g.dv_ER * (C*np.exp(-2.*g.dv_ER/VVT) - Ce) / (np.exp(-xxx)-1.) #: mole/meter**3/second  # dv_ER
        print("J1= ",J1)
        Jp     = (2*g.r*g.d_ER)/(g.N_A*dR2) * g.Omega_u * g.eta_p * C**2 / (C**2 + g.K_P**2) #: mole/meter**3/second # d_ER, N_A, Omega_u, eta_p, K_p
        print("Jp= ",Jp)
        #minf   =  I / (I + d_1)  : 1 # d_1
        #ninf   = C / (C + d_5) : 1 # d_5
        #hinf   =  d_2 * (I + d_1) / (d_2*(I +d_1) + (I+d_3)*C) : 1 # d_2, d_1, d_3

        # Must sum over two compartments. How to do this within this structure?

        rhs_C = Jr + J1 - Jp  #: mole / meter**3

        # Endoplasmic Reticulum Dynamnics
        rhs_Ce = Jp/(g.rho*g.Lambda) - (Jr + J1)  #: mole/meter**3

        #IP3 dynamics
        Jbeta  = 0  #*mole/meter**3/second  : mole/meter**3/second
        Jdelta = g.o_delta * (g.k_delta/(I+g.k_delta)) * (C**2/(C**2+g.K_delta**2)) #: mole/second  # not sure about units, o_delta, k_delta, K_delta
        J5P    = g.o_5P * (I/(I+g.K_5P)) #: mole/second # o_5P, K_5P
        J3K    = g.o_3K * (C**4/(C**4+g.K_D**4)) * (I/(I+g.K_3)) #: mole/second # o_3K, K_D, K_3

        # assume L constant (CHECK)
        Tot_I  = I                #: mole/meter**3  (sum of Ce in two half compartments in Brian2 code)
        coupling_I  = (4*g.D_I  / g.L**2) * (I  - I)  # : mole/second/meter**3 (summed) (CHECK)
        rhs_I = (Jbeta + Jdelta - J3K - J5P) / (g.Lambda*(1-g.rho)) #: mole/meter**3

        # Open Channel dynamics
        OmegaH = (g.Omega_2*(I+g.d_1) + g.O_2*(I+g.d_3)*C) / (I + g.d_3) #: Hz # Omega_2, O_2, d_1, d_3
        hinf   =  g.d_2 * (I + g.d_1) / (g.d_2*(I + g.d_1) + (I + g.d_3)*C) #: 1 # d_2, d_1, d_3
        rhs_h = OmegaH * (hinf - h) #: 1

        #-----------------------------
        # Handle diffusion terms, which involve terms from two adjacent compartments
        A1 = B1 = C1 = 0
        for i in [0,1]:
            A1 += (0.5 * g.F * s)/(g.Cm * g.V_T) * (dR2 / g.L)         #: meter**4 / mole (summed)
            # Default: fac=1. I think it should be 1000 or 0.001 because V_T is givin in mVolts
            fac = 1
            B1 += (1. - (fac * s * g.F * C[i])/(2. * g.Cm * g.V_T))            #: meter (summed)
            C1 += dR2 * C[i] / g.L                              #: mole / meter**2 (summed)

        print("A1,B1,C1= ", A1, B1, C1)

        nb_connections = 1
        CC0 = ((-B1 + np.sqrt(B1**2 + 4*C1*A1)) / (2.*A1)) / nb_connections #: mole / meter**3 (summed)
        V0 = (g.Vrest + (CC0 - g.Crest) * (g.F * s) / g.Cm) / nb_connections   #: volt  (summed)

        coupling_Ce = (4*g.D_CE / g.L**2) * (Ce[0]-Ce[1]) * np.array([1,-1])  #: mole/second/meter**3 (summed)  (CHECK) (most be Tot_Ce from each compartment
        electro_diffusion = -g.P_Ca * V / (g.R * g.V_T) * (Ce*np.exp(-2*V/g.V_T) - C) / (np.exp(-2*V/g.V_T) - 1) #: mole / second / meter**3
        coupling_electro = (4*g.D_C/g.L**2/g.V_T) * (CC0 + C) * (V0 - V)  #: mole/second/meter**3 (summed)
        coupling_C = (4*g.D_C / g.L**2) * (C[0] - C[1]) * np.array([1,-1]) #: mole/second/meter**3 (summed)  # ERROR. Cannot be zero.

        # All coupling terms
        rhs_C  += 1*coupling_C + 1.*coupling_electro + 0.*electro_diffusion #: mole / meter**3
        rhs_Ce += 1*coupling_Ce
        rhs_I  += 1*coupling_I   #: mole/meter**3

        print(rhs_C)
        print(rhs_Ce)
        print(rhs_I)
        print(rhs_h)

#----------------------------------------------------------------------
g = Parameters()

# Single compartment
rhs1 = RHS(g)
rhs2 = RHS(g)  # must clean this up. Probably compartment should only store state variables.
               # the RHS equations are really the same for all compartments, except diffusion terms.
               # not yet clear how to handle these.
C  = np.array([1,1])
Ce = np.array([1,1])
h  = np.array([1,1])
I  = np.array([1,1])
rhs1.rhs(C, Ce, h, I)  # arbitrary values of C, Ce, I, h
#----------------------------------------------------------------------

#Solver: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html
# If I advance wach rhs separately, as opposed to stacking them first, set h0 = hmin = hmax, to make sure
# that the stepsize used by the solver remains constant.
# I am not concerned with efficiency. This code is written to check against my Brian implementation.
