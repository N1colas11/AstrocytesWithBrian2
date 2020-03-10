using DifferentialEquations, Plots, BenchmarkTools
include("GlobalParameters.jl")
import .globalParameters_mod: globalParametersCorrect, globalParametersTest
getParams = globalParameters_mod.globalParametersEvan
#---------------------------------------------------------
# To run some parametric studies
# Add the name of a ODE solver as an argument, perhaps?
function run_problem(αs, βs)
    u0 = [1., 2.]
    tspan = (0., 3.)
    p = [1.0, 1.]
    αs = 1, 2
    βs = 3, 4
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

#run_problem(αs, βs)


function Hill(x, p, n)
    return @. x^n / (x^n + p^n)
end

lJp = []
lJr = []
lJ1 = []
lJleak = []
lJSERCA = []
lJIPR = []
lJβ = []
lCount = []
lTime = []

function Cresswell!(du, u, p, t)
    C  = @view(u[1:1:2])
    Ce = @view(u[3:1:4])
    I  = @view(u[5:1:6])
    h  = @view(u[7:1:8])
    push!(lCount, 1)
    #println("... Ce= ", Ce)

    Tot_C = C
    Tot_CE = Ce               #: mole/meter**3  (sum of Ce in two half compartments in Brian2 code)

    # C: calcium
    # Ce: calcium in ER
    # h : fraction of open channels
    # I : IP3 concentration

    # rhs_C, etc. Are arrays, created via implicit vectorization involving numpy arrays
    dR2 = p[:R]^2 - p[:r]^2  #: meter**2
    s = p[:R] + p[:r]
    V = @. p[:Vrest] + (p[:F]*s/p[:Cm]) * (C - p[:Crest]) #: volt
    VVT = -2. * V/p[:V_T] #: 1  # value of 0, which leads to a zero denominator
	#
    # Params:  r, P_r, p_open, dv_ER, P_CE, r, dv_ER, N_A, Omega_u, eta_p, K_P
    #Jr   = @. (2*p[:r]/dR2) * p[:P_r] * p[:p_open] * (Ce - C) #: mole/meter**3/second  # p_open, Pr
    #xxx  = @. -2. * p[:dv_ER]/VVT
    #J1   = @. (4 *p[:P_CE]/p[:r]) * VVT * p[:dv_ER] * (C * exp(-2 * p[:dv_ER] / VVT) - Ce) / (exp(-xxx) - 1.) #: mole/meter**3/second  # dv_ER
    #Jp   = @. (2 *p[:r]*p[:d_ER])/(p[:N_A]*dR2) * p[:Omega_u] * p[:eta_p] * C^2 / (C^2 + p[:K_P]^2) #: mole/meter**3/second # d_ER, N_A, Omega_u, eta_p, K_p


	minf    = @. Hill(I,p[:d1],1) * Hill(C,p[:d5],1)
    Jleak   = @. p[:vLeak] * (Ce-C)
    JIPR    = @. p[:vIPR] * (minf*h)^3 * (C-Ce)
    JSERCA  = @. p[:vSERCA] * Hill(C,p[:d5],2)

    # All three currents have the same units
    # (r/R^2) P_r p_open * C = (P_CE/r) * Volt * C  = (r * ER) (N_A * r^2) * Omega_u * eta_p
	#
    # push!(lJr, Jr)
    # push!(lJp, Jp)
    # push!(lJ1, J1)
	#
    push!(lJleak, Jleak)
    push!(lJIPR, JIPR)
    push!(lJSERCA, JSERCA)

	# Model used by Evan Cresswell for C, Ce
    #rhs_C  = @. JIPR - JSERCA + Jleak  #: mole / meter**3
    #rhs_Ce = @. -(1. / p[:rho_A]) * (JIPR - JSERCA + Jleak)  #: mole / meter**3
    # Endoplasmic Reticulum Dynamnics
    # J1 is electro coupling (turn off because of exponential)
    J_β = sin(2*π*t)^10  #*mole/meter**3/second  : mole/meter**3/second
    push!(lJβ, J_β)
	println("t= ", t, ",  typeof(t)= ", typeof(t))
	push!(lTime, t)

	J_δ = @. p[:o_δ]  * Hill(p[:k_δ], I, 1) * Hill(C, p[:K_δ], 2)
    J5P = @. p[:o_5P] * Hill(I, p[:K_5P], 1)  #: mole/second # o_5P, K_5P
    J3K = @. p[:o_3K] * Hill(C, p[:K_D], 4) * Hill(I, p[:K_3], 1) #: mole/second # o_3K, K_D, K_3

    # assume L constant (CHECK)
    #Tot_I  = I                #: mole/meter**3  (sum of Ce in two half compartments in Brian2 code)
    #coupling_I  = @. (4*p[:D_I]  / p[:L]^2) * (I  - I)  # : mole/second/meter**3 (summed) (CHECK)

    # Open Channel dynamics
    #OmegaH = @. (p[:Omega_2]*(I+p[:d_1]) + p[:O_2]*(I+p[:d_3])*C) / (I + p[:d_3]) #: Hz # Omega_2, O_2, d_1, d_3
    Q2     = @. p[:d2] * (I + p[:d1]) / (I + p[:d3])
	# hinf   = @. q2(p) / (Q2(p) + c) ;
	# tauh   = @. 1. / ( a2 * ( q2(p) + c  )  ) ;
	# OmegaH = 1/tauh
    OmegaH = @. p[:a2] * (Q2 + C)
    hinf   = @. p[:d_2] * (I + p[:d_1]) / (p[:d_2]*(I + p[:d_1]) + (I + p[:d_3])*C) #: 1 # d_2, d_1, d_3

    #-----------------------------
    # Handle diffusion terms, which involve terms from two adjacent compartments

    #print("A1,B1,C1= ", A1, B1, C1)

    nb_connections = 1

    dh  = @. OmegaH * (hinf - h) #: 1
    #dI  = @. (Jbeta + Jdelta - J3K - J5P) / (p[:Lambda]*(1-p[:rho])) #: mole/meter**3
    dI  = @. (J_β + J_δ - J3K - J5P)
    #println("dC= ", dC)

    #du .= dC, dCE, dI, dh  # .= to force update of list via reference
    # I need a more elegant approach
    # du[1:2] = dC
    # du[3:4] = dCE
    # du[5:6] = dI
    # du[7:8] = dh

    # simplified version for debuggin
    # Jr: Increase from .2 to 1 in time 1
    # J1: Increase from .2 to 150 in time 15
    # Jp: Constant solution
    # GE: divde Jp by Lambda (2020-03-09)
    dC  = @.                     (JIPR - JSERCA + Jleak)  #: mole / meter**3
    dCE = @. -(1. / p[:rho_A]) * (JIPR - JSERCA + Jleak)  #: mole / meter**3
    #  Jbeta: always zero
    #  Jdelta: goes to 3^4
    #  J3K: ``goes to -10^-12  (NOT GOOD)
    #  J5P: ``goes to -10^-13  (NOT GOOD) (goes to -.6 if I do not normalized by volume)

    println("Jbeta, Jdelta, J3K, J5P = ", J_β, J_δ, J3K, J5P)
    dI  = @. (1. *J_β + 1. * J_δ - 1. *J3K - 1. *J5P) # / (p[:Lambda]*(1-p[:rho])) #: mole/meter**3
    dh  = @. OmegaH * (hinf - h) #: 1
    #du .= dC, dCE, dI, dh  # .= to force update of list via reference
    #dC = [0,0]
    #dC = [0,0]
    du[1:2] = dC
    #du[1:2] = [0,0]
    du[3:4] = dCE #[0,0]
    #du[3:4] = [0,0]
    du[5:6] = dI # [0,0]  # blows up!
    du[7:8] = dh
end


nb_neurons = 2
eqs_per_neuron = 4
# Initial Conditions
function initialConditions()
    c0 = rand(nb_neurons)
    ce0 = rand(nb_neurons)
    h0 = rand(nb_neurons)
    I0 = rand(nb_neurons)
    # Initial condition vector
    u0 = vcat(c0,ce0,h0,I0)
end

function plotSol(sol)
	p1 = plot(sol.t, sol[1,:], label="C1")
	p1 = plot!(p1, sol.t, sol[2,:], label="C2")
	p2 = plot(sol.t, sol[3,:], label="Ce1")
	p2 = plot!(sol.t, sol[4,:], label="Ce2")
	p3 = plot(sol.t, sol[5,:], label="I1")
	p3 = plot!(p3,sol.t, sol[6,:], label="I2")
	p4 = plot(sol.t, sol[7,:], label="h1")
	p4 = plot!(p4,sol.t, sol[8,:], label="h2")
	plot(p1,p2,p3,p4, layout=(2,2))
end

u0 = initialConditions()
pars = getParams()
tspan = (0., 300.)  # make these real
#u0 = [.1, .2, .3, .4]``
# We have 8 equations. How to collect them?
prob = ODEProblem(Cresswell!, u0, tspan, pars)  # d is a parameter dictionary
@time sol = solve(prob, Tsit5())
#@time sol = solve(prob, Euler(), dt=1.e-7)
println(prob)
#interp = tspan[1]:(tspan[2]-tspan[1])/10.:tspan[2]
#sol_interp = sol(interp)
println("t= ", sol.t)
println("--------------------")
# Not using sol_interp correctly!
#plot(sol_interp[1,:])
println(u0)
plotSol(sol)
plot(lTime, lJβ[1:end-1], title="lJβ")
#plot!(p, sol.t, sol[7,:])
print(sol.t)
#print(lJ1[1,:])
print(size(lJ1[:,1]))
#plot(sol.t, lJ1)

#for i ∈ 1:8
    #print("i= ", i, ", sol= ", sol[i,:], "\n")
#end
lJβ
plot(time, lJβ)
