module globalParameters_mod
export globalParametersCorrect
export globalParametersEvan

#--------------------------------------------------------------
function globalParametersCorrect()
    d = Dict()
    d[:Rgas] = 8.31   # joule / kelvin / mole
    d[:D_I] = 5.3e-10 # Hz * meter^2
    d[:D_C] = 5.3e-10 # Hz * meter^2
    d[:D_CE] = 5.3e-10 # Hz * meter^2
    d[:Cm] = 1.e-6   # Farad * meter^-2
    d[:e] = 1.602176634e-19  # Amp * second
    d[:N_A] = 6.02214076e23 # mole^-1
    d[:k_B] = 1.38064852e-23  # meter^2 * Kg * Hz^2 / kelvin
    d[:F] = 96485.3329  # Amp * second / mole = Cb / mole
    d[:V_T] = 0.026  # voltage_cubic
    d[:P_Ca] = 4.46e-15  # meter * Hz
    d[:Crest] = 0.073e-6 # Molar
    d[:Vrest] = -80. # Volt
    d[:Ce] = 1.800  #  mole/meter**3]
    d[:Ce] = 1.e-4  #  mole/meter**3]  # WHICH ONE TO CHOOSE?

    d[:O_delta] = 0.6e-3  # Hz * mole/meter^3
    d[:k_delta] = 1.5e-3  # mole/meter^3
    d[:K_delta] = 0.1e-3  # mole/meter^3
    d[:K_D] = 0.7e-3  # mole/meter^3
    d[:K_3] = 1.0e-3  # mole/meter^3
    d[:O_3K] = 4.5e-3  # Hz*mole/meter^3
    d[:o_5P]   = 0.05   # Hz
    d[:K_5P] = 10.0e-3  # mole/meter^3
    d[:o_delta] = 0.6e-3  # mole/meter^3
    d[:k_delta] = 1.5e-3  # Hz * mole/meter^3
    d[:K_delta] = 0.1e-3   # mole/meter^3
    d[:Lambda] = 2100e-18  # meter^3
    d[:rho_A] = 0.18   # ER-to-cytoplasm volume ratio ?
    d[:rho] = d[:rho_A] / (1. + d[:rho_A])  # # ER-to-cytoplasm volume ratio ?

    # Multiply value of amplitudes by cytosolic volume (changed rho_A to rho)
    # It does not look as if I did the multiplication
    #self.o_delta = 0.6  # * umolar * Lambda * (1-rho) / second
    d[:o_delta] = 0.6e-3  # * mole/meter**3 * Lambda * (1-rho) / second
    d[:o_3K]    = 4.5e-3  # * mole/meter**3 * Lambda * (1-rho) / second

    d[:Omega_L] = 0.1   #/second         # Maximal rate of Ca^2+ leak from the ER
    d[:Omega_2] = 0.2   #/second      # IP_3R binding rate for Ca^2+ inhibition
    d[:Omega_u] = 0.2-3 #/second      # uptake/dissassociation constant

    # --- IP_3R kinectics
    d[:d_1] = 0.13e-3 #*mole/meter**3            # IP_3 binding affinity
    d[:d_2] = 1.05e-3 #*mole/meter**3            # Ca^2+ inactivation dissociation constant
    d[:O_2] = 0.2e3 #/(mole/meter**3)/second      # IP_3R binding rate for Ca^2+ inhibition
    d[:d_3] = 0.9434e-3 #*mole/meter**3          # IP_3 dissociation constant
    d[:d_5] = 0.08 #*mole/meter**3            # Ca^2+ activation dissociation constant

    d[:p_open] = 1
    d[:P_r]    = 1e-3 #*mole/meter**3/second
    d[:P_CE]    = 1e-3 #*mole/meter**3/second

    d[:eta_p] = 1
    d[:d_ER]  = 1e-6  #* meter
    d[:dv_ER] = 10e-3
    d[:K_P] = 0.05e-3

    d[:L] = 8e-6  # umeter
    d[:R] = 3e-6  # meter
    d[:r] = 2e-6  # meter
    d[:O_delta] = 0.6e-21  # mole/meter**3 * Lambda * (1-rho) / second
    d[:O_5P] = 0.05e-21   # mole/meter**3 * Lambda * (1-rho) / second
    d[:F] = 0.1  # ALREADY DEFINED?
    d[:D] = 0.05
    return d
end
#-----------------------------------------------------------
function globalParametersEvan()
    # Try to duplicate Evan's model as much as possible
    d = Dict()
    d[:Rgas] = 8.31   # joule / kelvin / mole
    d[:D_I] = 5.3e-10 # Hz * meter^2
    d[:D_C] = 5.3e-10 # Hz * meter^2
    d[:D_CE] = 5.3e-10 # Hz * meter^2
    d[:Cm] = 1.e-6   # Farad * meter^-2
    d[:e] = 1.602176634e-19  # Amp * second
    d[:N_A] = 6.02214076e23 # mole^-1
    d[:k_B] = 1.38064852e-23  # meter^2 * Kg * Hz^2 / kelvin
    d[:F] = 96485.3329  # Amp * second / mole = Cb / mole
    d[:V_T] = 0.026  # voltage_cubic
    d[:P_Ca] = 4.46e-15  # meter * Hz
    d[:Crest] = 0.073e-6 # Molar
    d[:Vrest] = -80. # Volt
    d[:Ce] = 1.800  #  mole/meter**3]
    d[:Ce] = 1.e-4  #  mole/meter**3]  # WHICH ONE TO CHOOSE?


    d[:O_delta] = 0.6e-3  # Hz * mole/meter^3
    d[:k_delta] = 1.5e-3  # mole/meter^3
    d[:K_delta] = 0.1e-3  # mole/meter^3
    d[:K_D] = 0.7e-3  # mole/meter^3
    d[:K_3] = 1.0e-3  # mole/meter^3
    d[:O_3K] = 4.5e-3  # Hz*mole/meter^3
    d[:o_5P]   = 0.05   # Hz
    d[:K_5P] = 10.0e-3  # mole/meter^3
    d[:o_delta] = 0.6e-3  # mole/meter^3
    d[:k_delta] = 1.5e-3  # Hz * mole/meter^3
    d[:K_delta] = 0.1e-3   # mole/meter^3
    d[:Lambda] = 2100e-18  # meter^3
    d[:rho_A] = 0.18   # ER-to-cytoplasm volume ratio ?
    d[:rho] = d[:rho_A] / (1. + d[:rho_A])  # # ER-to-cytoplasm volume ratio ?

    # Multiply value of amplitudes by cytosolic volume (changed rho_A to rho)
    # It does not look as if I did the multiplication
    #self.o_delta = 0.6  # * umolar * Lambda * (1-rho) / second
    d[:o_delta] = 0.6e-3  # * mole/meter**3 * Lambda * (1-rho) / second
    d[:o_3K]    = 4.5e-3  # * mole/meter**3 * Lambda * (1-rho) / second

    d[:Omega_L] = 0.1   #/second         # Maximal rate of Ca^2+ leak from the ER
    d[:Omega_2] = 0.2   #/second      # IP_3R binding rate for Ca^2+ inhibition
    d[:Omega_u] = 0.2-3 #/second      # uptake/dissassociation constant

    # --- IP_3R kinectics
    d[:d_1] = 0.13e-3 #*mole/meter**3            # IP_3 binding affinity
    d[:d_2] = 1.05e-3 #*mole/meter**3            # Ca^2+ inactivation dissociation constant
    d[:O_2] = 0.2e3 #/(mole/meter**3)/second      # IP_3R binding rate for Ca^2+ inhibition
    d[:d_3] = 0.9434e-3 #*mole/meter**3          # IP_3 dissociation constant
    d[:d_5] = 0.08 #*mole/meter**3            # Ca^2+ activation dissociation constant

    d[:p_open] = 1
    d[:P_r]    = 1e-3 #*mole/meter**3/second
    d[:P_CE]    = 1e-3 #*mole/meter**3/second

    d[:eta_p] = 1
    d[:d_ER]  = 1e-6  #* meter
    d[:dv_ER] = 10e-3
    d[:K_P] = 0.05e-3

    d[:L] = 8e-6  # umeter
    d[:R] = 3e-6  # meter
    d[:r] = 2e-6  # meter
    d[:O_delta] = 0.6e-21  # mole/meter**3 * Lambda * (1-rho) / second
    d[:O_5P] = 0.05e-21   # mole/meter**3 * Lambda * (1-rho) / second
    d[:F] = 0.1  # ALREADY DEFINED?
    d[:D] = 0.05

	# Evan parameters
	d[:vIPR] = 6  # sec^-1
	d[:vSERCA] = 4.4e-3  # mole m^-3 s^-1 (1 Molar = 1 Mole/liter = 10^3 mole/meter^3)
	d[:vLeak] = 0.11 # sec^-1
	d[:d1] = 0.13e-3 # mole m^-3
	d[:d5] = 0.08234e-3 # mole m^-3
	# End Evan Parameters
    return d
end
#---------------------------------------------
end
