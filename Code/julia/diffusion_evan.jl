using DifferentialEquations, Plots, BenchmarkTools
include("GlobalParameters.jl")
import .globalParameters_mod: globalParametersCorrect, globalParametersTest
include("global_functions.jl")
GF = globalFunctions_mod
GP = globalParameters_mod
getParams = GP.globalParametersEvan
#---------------------------------------------------------
# Global variables should be const for efficiency.
const nb_neurons = 2
const eqs_per_neuron = 5 # not used yet
const tot_nb_equations = nb_neurons * eqs_per_neuron # not used yet

function initialConditions()
    c0 = rand(nb_neurons)
    ce0 = rand(nb_neurons)
    h0 = rand(nb_neurons)
    I0 = rand(nb_neurons)
    prodβ0 = [0.,0.]
    u0 = vcat(c0,ce0,h0,I0, prodβ0)
end

function callback(u, t, integrator)
	C  = @view(u[1:1:2])  # Calcium
    Ce = @view(u[3:1:4])  # Calcium in ER
    I  = @view(u[5:1:6])  # IP3
    h  = @view(u[7:1:8])
	# How to I get my parameter array as a local variable?
	return currents = getCurrents(C, Ce, I, h, t, pars)
end

function Hill(x, p, n)
    return @. x^n / (x^n + p^n)
end

function getCurrents(C, Ce, I, h, t, p)
	minf    = @. Hill(I,p[:d1],1) * Hill(C,p[:d5],1)
	ninf    = @. Hill(C,p[:d5],1)
	hinf    = @. p[:d_2] * (I + p[:d_1]) / (p[:d_2]*(I + p[:d_1]) + (I + p[:d_3])*C) #: 1 # d_2, d_1, d_3
	Jleak   = @. p[:vLeak] * (Ce-C)
	JIPR    = @. p[:vIPR] * (minf*ninf*h)^3 * (C-Ce)
	JSERCA  = @. p[:vSERCA] * Hill(C,p[:d5],2)
	J_β     = sin(2*π*t)^10  #*mole/meter**3/second  : mole/meter**3/second
	J_δ     = @. p[:o_δ]  * Hill(p[:k_δ], I, 1) * Hill(C, p[:K_δ], 2)
	J5P     = @. p[:o_5P] * Hill(I, p[:K_5P], 1)  #: mole/second # o_5P, K_5P
	J3K     = @. p[:o_3K] * Hill(C, p[:K_D], 4) * Hill(I, p[:K_3], 1) #: mole/second # o_3K, K_D, K_3
	Q2      = @. p[:d2] * (I + p[:d1]) / (I + p[:d3])
	OmegaH  = @. p[:a2] * (Q2 + C)
	return J5P, J3K, J_β, J_δ, Jleak, JIPR, JSERCA, hinf, OmegaH, minf, ninf, Q2, t

end

function Cresswell!(du, u, p, t)
    C  = @view(u[1:1:2])  # Calcium
    Ce = @view(u[3:1:4])  # Calcium in ER
    I  = @view(u[5:1:6])  # IP3
    h  = @view(u[7:1:8])  # fraction of open channels
	prodβ = @view(u[9:1:10])  # Poisson pulses
    push!(lCount, 1)

    # rhs_C, etc. Are arrays, created via implicit vectorization involving numpy arrays
    dR2 = p[:R]^2 - p[:r]^2  #: meter**2
    s = p[:R] + p[:r]
    V = @. p[:Vrest] + (p[:F]*s/p[:Cm]) * (C - p[:Crest]) #: volt
    VVT = -2. * V/p[:V_T] #: 1  # value of 0, which leads to a zero denominator
	J5P, J3K, J_β, J_δ, Jleak, JIPR, JSERCA, hinf, OmegaH, minf, ninf, Q2, t_ = getCurrents(C, Ce, I, h, t, p)

	# Diffusion
	Diff_C  = [0.,0.]
	Diff_CE = [0.,0.]
	Diff_I  = @. p[:DI] * (I[2] - I[1]) * [1.,-1.]

	# dprodβ/dt = O_β * U[0,1; t_k] * ∑ δ(t-t_k)
	# Use callback to add the next Poisson pulse
	dprodβ = @. -p[:O_β] * prodβ * 200  # revisit the numerical coefficient
	J_β = prodβ

    # GE: divde Jp by Lambda (2020-03-09)
    dC  = @.                     (JIPR - JSERCA + Jleak)  + Diff_C #: mole / meter**3
    dCE = @. -(1. / p[:rho_A]) * (JIPR - JSERCA + Jleak)  + Diff_CE#: mole / meter**3
    dI  = @. (J_β + J_δ - J3K - J5P) + Diff_I
    dh  = @. OmegaH * (hinf - h) #: 1

    du[1:2]  = dC
    du[3:4]  = dCE #[0,0]
    du[5:6]  = dI # [0,0]  # blows up!
    du[7:8]  = dh
    du[9:10] = dprodβ
end

# Call back function to handle IP3 production
function condition(u,t,integrator)
	# return value is the condition that is >=0
	t - integrator.p[:tk]
end

function affect!(integrator)
	# In reality, update the time with Poisson interval
	#println("*** affect!: integrator.p[:tk] = ", integrator.p[:tk])
	integrator.p[:tk] = integrator.p[:tk] + 0.6
	#println("*** affect!: integrator.p[:tk] = ", integrator.p[:tk])
	# inefficient because memory allocation
	integrator.u[9:10] += @. integrator.p[:rhs_IP3_β] * rand() * [1., 1.]
end

function solve_problem()
	ccb = ContinuousCallback(condition, affect!)
	u0 = initialConditions()
	println("u0= ", u0)
	pars = getParams()
	tspan = (0., 5.)  # make these real
	#u0 = [.1, .2, .3, .4]``
	# We have 8 equations. How to collect them?
	prob = ODEProblem(Cresswell!, u0, tspan, pars)  # d is a parameter dictionary
	@time sol = solve(prob, Tsit5(), callback=ccb)
	#@time sol = solve(prob, Euler(), dt=1.e-7)
	println(prob)
	#interp = tspan[1]:(tspan[2]-tspan[1])/10.:tspan[2]
	#sol_interp = sol(interp)
	println("t= ", sol.u[10,:])  # glu spikes
	println("--------------------")
	# Not using sol_interp correctly!
	#plot(sol_interp[1,:])
	println(u0)
	GF.plotSol(sol)  # C, Ce, I, h, four plots
	plot(sol.t, sol[9,:]) # GLU pulses
	#plot!(p, sol.t, sol[7,:])
end

solve_problem()

plot(lTime, lJβ[1:end], title="lJβ")
