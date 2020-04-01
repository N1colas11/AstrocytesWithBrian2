using DifferentialEquations, Plots, BenchmarkTools
using Profile, ProfileView, TimerOutputs
include("GlobalParameters.jl")
import .globalParameters_mod: globalParametersCorrect, globalParametersTest
include("global_functions.jl")
const GF = globalFunctions_mod
const GP = globalParameters_mod
const getParams = GP.globalParametersEvan
#---------------------------------------------------------
# Global variables should be const for efficiency.
const nb_neurons = 2
const eqs_per_neuron = 5 # not used yet
const tot_nb_equations = nb_neurons * eqs_per_neuron # not used yet
const to = TimerOutput();

function initialConditions()
     #c0 = rand(nb_neurons)
    #ce0 = rand(nb_neurons)
     #h0 = rand(nb_neurons)
     #I0 = rand(nb_neurons)
     c0 = [.1, .1] .* 1.e-3
    ce0 = [.8, .8] .* 1.e-3
     #h0 = rand(nb_neurons)
     h0 = [0.1, 0.1]
     I0 = [0.1, 0.1] * 1.e-1
     #I0 = rand(nb_neurons)
    prodβ0 = [0.,0.]
    u0 = vcat(c0,ce0,h0,I0, prodβ0)
end


#@inline function Hill(x, p, n)
    #return @. x^n / (x^n + p^n)
#end
@inline Hill(x, p, n) = @. x^n / (x^n + p^n)

function getCurrents(C, Ce, I, h, t, p)
  @timeit to "getCurrents" begin
    minf    = @. Hill(I,p[:d1],1)
	ninf    = @. Hill(C,p[:d5],1)
	hinf    = @. p[:d_2] * (I + p[:d_1]) / ((I + p[:d_1]) + (I + p[:d_3])*C) #: 1 # d_2, d_1, d_3
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
end

function Cresswell!(du, u, p, t)
    @timeit to "Cresswell!" begin # Must be at least one space or won't work
    C  = @view(u[1:1:2])  # Calcium
    Ce = @view(u[3:1:4])  # Calcium in ER
    I  = @view(u[5:1:6])  # IP3
    h  = @view(u[7:1:8])  # fraction of open channels
	prodβ = @view(u[9:1:10])  # Poisson pulses

	# CURRENTS
    minf    = @. Hill(I,p[:d1],1)
	ninf    = @. Hill(C,p[:d5],1)
	hinf    = @. p[:d_2] * (I + p[:d_1]) / ((I + p[:d_1]) + (I + p[:d_3])*C) #: 1 # d_2, d_1, d_3
	Jleak   = @. p[:vLeak] * (Ce-C)
	JIPR    = @. p[:vIPR] * (minf*ninf*h)^3 * (C-Ce)
	JSERCA  = @. p[:vSERCA] * Hill(C,p[:d5],2)
	J_β     = sin(2*π*t)^10  #*mole/meter**3/second  : mole/meter**3/second
	J_δ     = @. p[:o_δ]  * Hill(p[:k_δ], I, 1) * Hill(C, p[:K_δ], 2)
	J5P     = @. p[:o_5P] * Hill(I, p[:K_5P], 1)  #: mole/second # o_5P, K_5P
	J3K     = @. p[:o_3K] * Hill(C, p[:K_D], 4) * Hill(I, p[:K_3], 1) #: mole/second # o_3K, K_D, K_3
	Q2      = @. p[:d2] * (I + p[:d1]) / (I + p[:d3])
	OmegaH  = @. p[:a2] * (Q2 + C)


    # rhs_C, etc. Are arrays, created via implicit vectorization involving numpy arrays
    dR2 = p[:R]^2 - p[:r]^2  #: meter**2
    s = p[:R] + p[:r]
    V = @. p[:Vrest] + (p[:F]*s/p[:Cm]) * (C - p[:Crest]) #: volt
    VVT = -2. * V/p[:V_T] #: 1  # value of 0, which leads to a zero denominator
	#J5P, J3K, J_β, J_δ, Jleak, JIPR, JSERCA, hinf, OmegaH, minf, ninf, Q2, t_ = getCurrents(C, Ce, I, h, t, p)

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
    #println("J_β= ", J_β, ",  J_δ= ", J_δ, ", J3K= ", J3K, ", J5P= ", J5P)
    dh  = @. OmegaH * (hinf - h) #: 1

    du[1:2]  = dC
    du[3:4]  = dCE #[0,0]
    du[5:6]  = dI # [0,0]  # blows up!
    du[7:8]  = dh
    du[9:10] = dprodβ
end
end

function solve_problem(final_time::Float64=1.)
	# Call back function to handle IP3 production
	function condition(u,t,integrator)
		# return value is the condition that is >=0
		integrator.p[:tk] - t
	end

	function affect!(integrator)
		# In reality, update the time with Poisson interval
		integrator.p[:tk] = popfirst!(integrator.p[:events_iterator])
		#print("affect!(integrator), tk:  ", integrator.p[:tk])
		# inefficient because memory allocation
		integrator.u[9:10] += @. integrator.p[:rhs_IP3_β] * rand() * [1., 1.]
	end

	ccb = ContinuousCallback(condition, affect!)
	u0 = initialConditions()
	#println("u0= ", u0)
	pars = getParams()
	tspan = (0., final_time)  # make these real
	λ = 0.3
	tks = GF.generatePoissonEvents(2000, λ)
	#println("Poisson events: ", tks[1:10], "...")
	pars[:events] = tks
	pars[:events_iterator] = Iterators.Stateful(tks)
	pars[:tk] = popfirst!(pars[:events_iterator])
	#println("tk= ", pars[:tk])
	#println("events_iterator= ", pars[:events_iterator])
	#u0 = [.1, .2, .3, .4]``
	# We have 8 equations. How to collect them?
	#println("ODEProblem")
	prob = ODEProblem(Cresswell!, u0, tspan, pars)  # d is a parameter dictionary
	#println("solve")
	@timeit to "solve" @time sol = solve(prob, Tsit5(), callback=ccb)
	#@time sol = solve(prob, Euler(), dt=1.e-7)
	#println("GF.plotSol")
	#plot = GF.plotSol(sol)  # C, Ce, I, h, four plots
	#display(plot)
	#plot(sol.t, sol[9,:]) # GLU pulses
	#plot!(p, sol.t, sol[7,:])
	return sol
end

reset_timer!(to) # of TimerOutputs
sol = solve_problem(10.); println(to)

plot(sol.t, sol[9,:])
GF.plotSol(sol)
