using DifferentialEquations, Plots, BenchmarkTools
using Profile, ProfileView, TimerOutputs
using DataFrames
#include("GlobalParameters.jl")
#import .globalParameters_mod: globalParametersCorrect, globalParametersTest
#include("global_functions.jl")
#const GF = globalFunctions_mod
#const GP = globalParameters_mod
#const getParams = GP.globalParametersEvan

function original_model()
	# Taken from R code. This is what Nicolas did, but I want my own interpretation


end

function globalParameters()
  d = Dict{Symbol,Float64}()
  d[:base_mortality_ratio] = 1 / 30.
  d[:age_mortality_ratio]  = 1 / 50.

  d[:step_up] = 180
  d[:overcrowd_cost] = 2
  d[:overcrowd_thresh] = 0.5e6
  d[:n_step] = 1000
  d[:step_size] = 1
  d[:plot] = 0  # equivalent to false

  d[:year]         = 365    # year, time point at which transmission levels return to normal values
  d[:revertR]      = 2.8    # COVID-19 R0 (the average number of individuals infected by a single individual in a completely susceptible population)"
  d[:recovery_time] = 14     # recovery time in days. Taken to be 14 day
  d[:hosp_mortality_ratio]      = 10     # the hospitalization/mortality ratio

  d[:N]            = 300e6  # Population size  (300 million in the US)
  d[:Io]           = 1000   # Population Infected over 65
  d[:Iu]           = 1000   # Population Infected under 65
  d[:ρo]           = 0.15   # fraction of population over 65
  d[:Io]           = 1000   # number of infected individuals
  d[:Iy]           = 1000
  d[:αu]           = .1
  d[:αo]           = .1
  d[:βuu]          = .1     # under 65 / under 65 infection rate
  d[:γ0]           = .1     # over 65 recovery rate
  d[:δ0]           = .1     # over 65 mortality rate
  d[:Sf]           = 1.0    # fraction of total population suceptible to infection
  return d
end

function derivedParameters(d)
	d[:ρu]  = 1. - d[:ρo]  # fraction of population under 65
	d[:γu]  = d[:γo] = d[:γ0]
	d[:δu]  = d[:δo] = d[:δ0]
	d[:β]   = d[:βuo] = d[:βou] = d[:βoo] = d[:βuu]     # under 65 / over 65 infection rate
	#d[:αu]   = d[:γ0] * d[:ρ0] # Recovery rate + Death rate for old people
	#ay=gy+Dy  #Recovery rate + Death rate for young people
	d[:α0]   = d[:αu]*d[:ρu] + d[:αo]*d[:ρo]
	d[:go] = d[:gu] = 1. / d[:recovery_time] #recovery rate
    d[:No0] = d[:N] * d[:ρo]  # initial population over 65   (called No0 and Ny0)
    d[:Nu0] = d[:N] * d[:ρu]  # initial population under 65

	# Initial Conditions. But I will call the initial condition function later,
	# which will use these paramters
	d[:Mo0] = 0                  # number of dead at t=0
	d[:Mu0] = 0
	d[:So0] = d[:No0] * d[:Sf]   # Total number of susceptibles at t=0
	d[:Su0] = d[:Nu0] * d[:Sf]
	d[:Ro0] = d[:No0] - d[:So0]  # Total number of recovered at =t0
	d[:Ru0] = d[:Nu0] - d[:Su0]
	# Initial fraction recovered population
	d[:Io0] = d[:No0] - d[:So0] - d[:Ro0] - d[:Mo0]
	d[:Iu0] = d[:Nu0] - d[:Su0] - d[:Ru0] - d[:Mu0]
	println("sum (I+S+R+M) old+young: ",
	   d[:Io0]+d[:Iu0]+d[:Ro0]+d[:Ru0]+d[:So0]+d[:Su0]+d[:Mo0]+d[:Mu0])
	println("gordon: ",
	    "Frances")
end

function initialConditions(p::Dict{Symbol,Float64})
	# args p are parameters
     I0 = p[:ρo]
     S0 = p[:ρo] * p[:N]
     R0 = .1
     M0 = .1
     i0 = p[:ρu]
     s0 = p[:ρu] * p[:N]
     r0 = .1
     m0 = .1
     uinit = vcat(I0,S0,R0,M0,i0,s0,r0,m0)
     return uinit
end


function Covid!(du, u, p, t)
    #@timeit to "Covid!" begin # Must be at least one space or won't work
    I  = @view(u[1])[1]  #
    S  = @view(u[2])[1]  #
    R  = @view(u[3])[1]  #
    M  = @view(u[4])[1]  #
    i  = @view(u[5])[1]  #
    s  = @view(u[6])[1]  #
    r  = @view(u[7])[1]  #
    m  = @view(u[8])[1]  #
    #println("I,S,i,s= ", I,S,i,s)
    #println("I,S,i,s= $I, $S, $i, $s")

	  αu = p[:γu] + p[:δu]
	  αo = p[:γo] + p[:δo]

	  du[1] =  p[:βuo]*S*i/p[:N] + p[:βoo]*S*I/p[:N] - p[:αo]*I   # I
	  du[2] = -p[:βuo]*S*i/p[:N] - p[:βoo]*S*I/p[:N]              # S
	  du[3] = +p[:γo]*I                                           # R
	  du[4] = +p[:δo]*I                                           # M
	  du[5] =  p[:βuu]*s*i/p[:N] + p[:βuo]*s*I/p[:N] - p[:αu]*I   # i
	  du[6] = -p[:βuu]*s*i/p[:N] - p[:βuo]*s*I/p[:N]              # s
	  du[7] = +p[:γu]*I                                           # r
	  du[8] = +p[:δu]*I                                           # m
	  old   = I+S+R+M
	  young = i+s+r+m
	  println("young+old: ", young+old)

# d/dt [I + S + R + M] = (-p[:αo]*I + p[:γo] + p[:δo]) * I = 0
# d/dt [i + s + r + m] = (-p[:αu]*I + p[:γu] + p[:δu]) * I = 0
#
# Old people stay old; young people stay young, so the balance
#   above holds for each class individually.
# Therefore, αo = γo + δo
#    and     αu = γu + δu
#
#    I + S + R + M = const = p[:ρo]
#    i + s + r + m = const = p[:ρu]

end

pp = globalParameters()
derivedParameters(pp)
u0 = initialConditions(pp)
tspan = (0., 365.)
#tspan = (0., 5.)


prob = ODEProblem(Covid!, u0, tspan, pp)  # d is a parameter dictionary
@time sol = solve(prob, Tsit5())
#@time sol = solve(prob, Euler(), dt=.1)
labels=[:I,:S,:R,:M,:i,:s,:r,:m]
p = plot()
for (i,label) ∈ enumerate(labels)
	#print(i)
	plot!(p, sol.t, sol[i,:], label=label)
end
#print(sol)
p
