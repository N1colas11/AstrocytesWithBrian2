using DifferentialEquations, Plots, BenchmarkTools


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

function globalParameters()
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


function CIh!(du, u, p, t)
    #C, Ce, I, h = u
    C = @view(u[1:1:2])
    Ce = @view(u[3:1:4])
    I = @view(u[5:1:6])
    h = @view(u[7:1:8])
    println("... Ce= ", Ce)

    Tot_C = C
    Tot_CE = Ce               #: mole/meter**3  (sum of Ce in two half compartments in Brian2 code)

    # C: calcium
    # Ce: calcium in ER
    # h : fraction of open channels
    # h : IP3 concentration

    # rhs_C, etc. Are arrays, created via implicit vectorization involving numpy arrays
    dR2 = p[:R]^2 - p[:r]^2  #: meter**2
    s = p[:R] + p[:r]
    V = p[:Vrest] .+ (p[:F]*s/p[:Cm]) * (C .- p[:Crest]) #: volt
    VVT = -2. * V/p[:V_T] #: 1  # value of 0, which leads to a zero denominator

    J_C = diff_C   = 0
    J_Ce = diff_Ce = 0
    J_I = diff_I   = 0

    Jr     = @. (2*p[:r]/dR2) * p[:P_r] * p[:p_open] * (Ce - C) #: mole/meter**3/second  # p_open, Pr
    xxx   = @. -2. * p[:dv_ER]/VVT
    println("xxx= ", xxx)
    println("VVT= ", VVT)
    println("C= ", C)
    println("Ce= ", Ce)
    # J1 has 4 elements instead of 2!! Not possible
    J1   = @. (4 *p[:P_CE]/p[:r]) * VVT * p[:dv_ER] * (C * exp(-2 * p[:dv_ER] / VVT) - Ce) / (exp(-xxx) - 1.) #: mole/meter**3/second  # dv_ER
    Jp   = @. (2 *p[:r]*p[:d_ER])/(p[:N_A]*dR2) * p[:Omega_u] * p[:eta_p] * C*2 / (C^2 + p[:K_P]^2) #: mole/meter**3/second # d_ER, N_A, Omega_u, eta_p, K_p

    println("Jr= ", Jr)
    println("J1= ", J1)
    println("Jp= ", Jp)
    rhs_C = @. Jr + J1 - Jp  #: mole / meter**3

    # Endoplasmic Reticulum Dynamnics
    rhs_Ce = @. Jp/(p[:rho]*p[:Lambda]) - (Jr + J1)  #:
    Jbeta  = 0  #*mole/meter**3/second  : mole/meter**3/second
    Jdelta = @. p[:o_delta] * (p[:k_delta]/(I+p[:k_delta])) * (C^2/(C^2+p[:K_delta]^2)) #: mole/second  # not sure about units, o_delta, k_delta, K_delta
    J5P    = @. p[:o_5P] * (I/(I+p[:K_5P])) #: mole/second # o_5P, K_5P
    J3K    = @. p[:o_3K] * (C^4/(C^4+p[:K_D]^4)) * (I/(I+p[:K_3])) #: mole/second # o_3K, K_D, K_3

    # assume L constant (CHECK)
    Tot_I  = I                #: mole/meter**3  (sum of Ce in two half compartments in Brian2 code)
    coupling_I  = @. (4*p[:D_I]  / p[:L]^2) * (I  - I)  # : mole/second/meter**3 (summed) (CHECK)

    # Open Channel dynamics
    OmegaH = @. (p[:Omega_2]*(I+p[:d_1]) + p[:O_2]*(I+p[:d_3])*C) / (I + p[:d_3]) #: Hz # Omega_2, O_2, d_1, d_3
    hinf   = @. p[:d_2] * (I + p[:d_1]) / (p[:d_2]*(I + p[:d_1]) + (I + p[:d_3])*C) #: 1 # d_2, d_1, d_3

    #-----------------------------
    # Handle diffusion terms, which involve terms from two adjacent compartments
    A1 = 0
    B1 = 0
    C1 = 0

    for i ∈ [1,2]
        A1 += @. (0.5 * p[:F] * s)/(p[:Cm] * p[:V_T]) * (dR2 / p[:L])         #: meter**4 / mole (summed)
        # Default: fac=1. I think it should be 1000 or 0.001 because V_T is givin in mVolts
        fac = 1
        println(C[i])
        B1 += @. (1. - (fac * s * p[:F] * C[i])/(2. * p[:Cm] * p[:V_T]))            #: meter (summed)
        C1 += @. dR2 * C[i] / p[:L]                              #: mole / meter**2 (summed)
    end

    print("A1,B1,C1= ", A1, B1, C1)

    nb_connections = 1
    CC0               = @. ((-B1 + sqrt(B1^2 + 4*C1*A1)) / (2*A1)) / nb_connections #: mole / meter**3 (summed)
    V0                = @. (p[:Vrest] + (CC0 - p[:Crest]) * (p[:F] * s) / p[:Cm]) / nb_connections   #: volt  (summed)
    coupling_Ce       = @. (4*p[:D_CE] / p[:L]^2) * (Ce[1]-Ce[2]) * (1,-1)   #: mole/second/meter**3 (summed)  (CHECK) (most be Tot_Ce from each compartment
    electro_diffusion = @. -p[:P_Ca] * V / (p[:R] * p[:V_T]) * (Ce*exp(-2*V/p[:V_T]) - C) / (exp(-2*V/p[:V_T]) - 1) #: mole / second / meter**3
    coupling_electro  = @. (4*p[:D_C]/p[:L]^2/p[:V_T]) * (CC0 + C) * (V0 - V)  #: mole/second/meter**3 (summed)
    coupling_C        = @. (4*p[:D_C] / p[:L]^2) * (C[1] - C[2]) * (1,-1) #: mole/second/meter**3 (summed)  # ERROR. Cannot be zero.

    dh  = @. OmegaH * (hinf - h) #: 1
    dI  = @. (Jbeta + Jdelta - J3K - J5P) / (p[:Lambda]*(1-p[:rho])) #: mole/meter**3
    dC  = @. 1*coupling_C + 1*coupling_electro + 0*electro_diffusion + Jr + J1 - Jp  #: mole / meter**3
    dCE = @. 1*coupling_Ce + Jp/(p[:rho]*p[:Lambda]) - (Jr + J1)  #: mole/meter**3
    print("dC= ", dC)

    #du .= dC, dCE, dI, dh  # .= to force update of list via reference
    # I need a more elegant approach
    du[1:2] = dC
    du[3:4] = dCE
    du[5:6] = dI
    du[7:8] = dh
end



nb_neurons = 2
eqs_per_neuron = 4
# Initial Conditions
c0 = rand(nb_neurons)
ce0 = rand(nb_neurons)
h0 = rand(nb_neurons)
I0 = rand(nb_neurons)
u0 = vcat(c0,ce0,h0,I0)

pars = globalParameters()
tspan = (0., 1.)  # make these real
#u0 = [.1, .2, .3, .4]
# We have 8 equations. How to collect them?
#CIh!(du, u0, p, t)
prob = ODEProblem(CIh!, u0, tspan, pars)  # d is a parameter dictionary
println(prob)
sol = solve(prob, Tsit5())
println(prob)
println("--------------------")

u0 = [1., 2.]
tspan = (0., 3.)
p = [1.0, 1.]
αs = 1, 2
βs = 3, 4

function solve_CIh!(u0, tspan, d)
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
