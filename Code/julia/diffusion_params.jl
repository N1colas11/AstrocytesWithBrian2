using DifferentialEquations
using Plots
using BenchmarkTools

function coupled!(du, u, p, t)
    α, β, D1, D2 = p
    c1, c2 = u
    dc1  = α*c1 + D1 * (c2-c1)
    dc2 = β*c2  + D2 * (c1-c2)
    du .= [dc1, dc2]  # note the broadcast operator!
end

function Hill(x, p, n)
    return @. x^n / (x^n + p^n)
end

x = rand(10000)

@btime a = Hill(x, 1., 2)



function CIh(du, u, p, t)
    Tot_C = C
    Tot_CE = Ce               #: mole/meter**3  (sum of Ce in two half compartments in Brian2 code)

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

    J_C = diff_C   = 0
    J_Ce = diff_Ce = 0
    J_I = diff_I   = 0

    rhs_C  = J_C + diff_C
    rhs_Ce = J_Ce + diff_Ce
    rhs_I  = J_I + diff_I

    xxx   = -2.*d[:dv_ER]/d[:VVT]
    J1    = (4*p[:P_CE]/p[:r]*p[:VVT]) * p[:dv_ER] * (C*np.exp(-2.*p[:dv_ER]/p[:VVT]) - Ce) / (np.exp(-xxx)-1.) #: mole/meter**3/second  # dv_ER
    print("J1= ",J1)
    Jp    = (2*g.r*g.d_ER)/(g.N_A*dR2) * g.Omega_u * g.eta_p * C**2 / (C**2 + g.K_P**2) #: mole/meter**3/second # d_ER, N_A, Omega_u, eta_p, K_p

    rhs_C = Jr + J1 - Jp  #: mole / meter**3

    # Endoplasmic Reticulum Dynamnics
    rhs_Ce = Jp/(g.rho*g.Lambda) - (Jr + J1)  #:
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

    dI = (Jbeta + Jdelta - J3K - J5P) / (Lambda*(1-rho)) + 1*coupling_CE   : mole/meter**3
    dC = 1*coupling_C + 1.*coupling_electro + 0.*electro_diffusion + Jr + J1 - Jp  : mole / meter**3
    dCE = 1*coupling_CE + Jp/(rho*Lambda) - (Jr + J1)  : mole/meter**3


    du .= dC, dCE, dI, dh
end

function globalParameters()
    d = Dict()
    d[:Rgas] = 8.31   # joule / kelvin / mole
    d[:D_I] = 5.3e-10 # Hz * meter^2
    d[:D_C] = 5.3e-10 # Hz * meter^2
    println("D_C: ", d[:D_C])
    d[:Cm] = 1.e-6   # Farad * meter^-2
    d[:e] = 1.602176634e-19  # Amp * second
    d[:N_A] = 6.02214076e23 # mole^-1
    d[:k_B] = 1.38064852e-23  # meter^2 * Kg * Hz^2 / kelvin
    d[:F] = 96485.3329  # Amp * second / mole = Cb / mole
    d[:V_T] = 0.026  # voltage_cubic
    d[:P_Ca] = 4.46e-15  # meter * Hz
    d[:Crest] = 0.073e-6 # Molar
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

d = globalParameters()
println(d)


u0 = [1., 2.]
tspan = (0., 3.)
p = [1.0, 1.]
αs = 1, 2
βs = 3, 4

function solve_CIh!(u0, tspan, d)
    prob = ODEProblem(CIh!, u0, tspan, d)  # d is a parameter dictionary
    solve(prob, Tsit5())
end



function run_problem(αs, βs)
    # Run all α, β pairs
    println(αs)
    sols = Any[]
    #map(αs,βs) do α, β
    for α ∈ αs
        for β ∈ βs
            p1 = [α, β, 8., 5.]  # α, β, D1, D2
            prob = ODEProblem(coupled!, u0, tspan, p1)
            sol = solve(prob, Tsit5())
            push!(sols, sol)
        end
    end
    p = plot(sols[1])
    for i ∈ 2:4
        plot!(p, sols[i])   # first argument: update plot
    end
    #display(p)   # necessary if plot is in a function
    display(p)    # or just p if it is the action of this method
end

run_problem(αs, βs)

t = [-1,0,1]
a = test_plot()
plot(t, a[1])
plot!(t, a[2])

plot(a[2])
