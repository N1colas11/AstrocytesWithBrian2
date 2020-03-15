using DifferentialEquations, Plots, BenchmarkTools
using DataFrames, Pandas
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


#cols = [:minf, :Jleak, :JIPR, :JSERCA, :J_β, :J_δ, :J5P, :J3K, :Q2, :OmegaH, :hinf]
#df = DF.DataFrame(VF, cols)
#let df = DF.DataFrame(VF, cols)
	function getCurrents(C, Ce, I, h, t, p)
		minf    = @. Hill(I,p[:d1],1) * Hill(C,p[:d5],1)
		Jleak   = @. p[:vLeak] * (Ce-C)
		JIPR    = @. p[:vIPR] * (cu[:minf]*h)^3 * (C-Ce)
		JSERCA  = @. p[:vSERCA] * Hill(C,p[:d5],2)
		J_β     = [sin(2*π*t)^10, cos(2*π*t)]  #*mole/meter**3/second  : mole/meter**3/second  (just to have a vector)
		J_δ     = @. p[:o_δ]  * Hill(p[:k_δ], I, 1) * Hill(C, p[:K_δ], 2)
		J5P     = @. p[:o_5P] * Hill(I, p[:K_5P], 1)  #: mole/second # o_5P, K_5P
		J3K     = @. p[:o_3K] * Hill(C, p[:K_D], 4) * Hill(I, p[:K_3], 1) #: mole/second # o_3K, K_D, K_3
		Q2      = @. p[:d2] * (I + p[:d1]) / (I + p[:d3])
		OmegaH  = @. p[:a2] * (cu[:Q2] + C)
		hinf    = @. p[:d_2] * (I + p[:d_1]) / (p[:d_2]*(I + p[:d_1]) + (I + p[:d_3])*C) #: 1 # d_2, d_1, d_3
		cols = (minf, Jleak, JIPR, JSERCA, J_β, J_δ, J5P, J3K, Q2, OmegaH, hinf)
		return cols
		# return a named tuple
		#V = Vector{Any}
		#named_tuple = NamedTuple(cols, Tuple(V,V,V,V,V,V,V,V,V,V,V))
		#return cols, ()

		#if length(df) == 0
            #VF = Vector{Any}
            #eltype = [VF, VF, VF, VF, VF, VF, VF, VF, VF, VF, VF]
			#df = DataFrames.DataFrame(cols)
			#print(typeof(df))
		#else
		    #push!(df, cols)
		#return df
	end
end

function Cresswell!(du, u, p, t)
    C  = @view(u[1:1:2])  # Calcium
    Ce = @view(u[3:1:4])  # Calcium in ER
    I  = @view(u[5:1:6])  # IP3
    h  = @view(u[7:1:8])  # fraction of open channels
    push!(lCount, 1)

    # rhs_C, etc. Are arrays, created via implicit vectorization involving numpy arrays
    dR2 = p[:R]^2 - p[:r]^2  #: meter**2
    s = p[:R] + p[:r]
    V = @. p[:Vrest] + (p[:F]*s/p[:Cm]) * (C - p[:Crest]) #: volt
    VVT = -2. * V/p[:V_T] #: 1  # value of 0, which leads to a zero denominator

	J5P, J3K, J_β, J_δ, Jleak, JIPR, JSERCA, hinf, OmegaH, minf, Q2 = getCurrents(C, Ce, I, h, t, p)
	#cu = getCurrents(C, Ce, I, h, t, p)

    # GE: divde Jp by Lambda (2020-03-09)
    #dC  = @.                     (cu[:JIPR] - cu[:JSERCA] + cu[:Jleak])  #: mole / meter**3
    #dCE = @. -(1. / p[:rho_A]) * (cu[:JIPR] - cu[:JSERCA] + cu[:Jleak])  #: mole / meter**3
    #dI  = @. (cu[:J_β] + cu[:J_δ] - cu[:J3K] - cu[:J5P])
    #dh  = @. cu[:OmegaH] * (cu[:hinf] - h) #: 1
    dC  = @.                     (JIPR - JSERCA + Jleak)  #: mole / meter**3
    dCE = @. -(1. / p[:rho_A]) * (JIPR - JSERCA + Jleak)  #: mole / meter**3
    dI  = @. (J_β + J_δ - J3K - J5P)
    dh  = @. OmegaH * (hinf - h) #: 1

    du[1:2] = dC
    du[3:4] = dCE #[0,0]
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

function callback(u, t, integrator)
	C  = @view(u[1:1:2])  # Calcium
    Ce = @view(u[3:1:4])  # Calcium in ER
    I  = @view(u[5:1:6])  # IP3
    h  = @view(u[7:1:8])
	println("iside callback")
	# How to I get my parameter array as a local variable?
	#J5P, J3K, J_β, J_δ, Jleak, JIPR, JSERCA, hinf, OmegaH, minf, Q2 = getCurrents(C, Ce, I, h, t, p)
	currents = getCurrents(C, Ce, I, h, t, pars)
	println("currents= ", currents)
	try
		println("***** inside try ****")
	    push!(df, currents)
		println("try, df= ", df)
    catch
		V = Vector{Any}
		S = Float64
		println("**** inside catch *****")
        #df = DF.DataFrame([V,V,V,V,V,V,V,V,V,V,V], [:J_5P, :J_3K, :J_β, :J_δ, :Jleak, :JIPR, :JSERCA, :hinf, :OmegaH, :minf, :Q2])
		V = Vector{Any}
		cols = (:minf, :Jleak, :JIPR, :JSERCA, :J_β, :J_δ, :J5P, :J3K, :Q2, :OmegaH, :hinf)
		#tuple1 = Tuple{V,V,V}
		#tuple = Tuple{V,V,V,V,V,V,V,V,V,V,V}
		#named_tuple = NamedTuple(cols, Tuple{V,V,V,V,V,V,V,V,V,V,V})
		named_tuple = NamedTuple{cols, Tuple{V,V,V,V,V,V,V,V,V,V,V}}
		push!(DF.DataFrame(), named_tuple)
		println("catch, df= ", df)
	    println("catch, currents= ", currents)
	    println("catch, after push!")
		println("df.names= ", names(df))
		println("catch, after push, df= ", df)
		push!(df, currents)
		println("catch, df= ", df)
		println("catch, df= ", df.J_5P)
	end
	println("after try/catch")
	#return currents = getCurrents(C, Ce, I, h, t, pars)
	return df
end

a = [3,4,5]; println(typeof(a))
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

saved_values = SavedValues(Float64, Vector{Any})
#saved_values = SavedValues(Float64, Dict{Any,Any})
cb = SavingCallback(callback, saved_values)

u0 = initialConditions()
pars = getParams()
tspan = (0., 10.)  # make these real
#u0 = [.1, .2, .3, .4]``
# We have 8 equations. How to collect them?
prob = ODEProblem(Cresswell!, u0, tspan, pars)  # d is a parameter dictionary
@time sol = solve(prob, Tsit5(), callback=cb)
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
plot(lTime, lJβ[1:end], title="lJβ")
#plot!(p, sol.t, sol[7,:])
print(sol.t)

function tst()
	V = Vector{Any}
	cols = (:minf, :Jleak, :JIPR, :JSERCA, :J_β, :J_δ, :J5P, :J3K, :Q2, :OmegaH, :hinf)
	vals = ([V for i ∈ 1:7])
	tuple = Tuple([V for i ∈ 1:11])
	tuple = Tuple([[0.,0.] for i ∈ 1:11])
	println(tuple)
	df = DF.DataFrame(tuple)
	append!(df, tuple)
	println(df)
	#named_tuple = NamedTuple{cols, Tuple{V,V,V,V,V,V,V,V,V,V,V}}
	#push!(DF.DataFrame(), named_tuple)
end
tst()
