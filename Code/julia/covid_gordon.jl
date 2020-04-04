using DifferentialEquations, Plots, BenchmarkTools,
     Profile, ProfileView, TimerOutputs,
     DataFrames


function globalParameters()
  d = Dict{Symbol,Float64}()
  d[:mult_init] =
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
  d[:βuu]          = .1     # under 65 / under 65 infection rate
  d[:γ0]           = .1     # over 65 recovery rate (assumed constant)
  d[:δ0]           = .1     # over 65 mortality rate (assumed constant)
  d[:Sf]           = 1.0    # fraction of total population suceptible to infection
  return d
end

function derivedParameters(d)
	d[:ρu]  = 1. - d[:ρo]  # fraction of population under 65
	d[:γu]  = d[:γo] = d[:γ0]
	d[:δu]  = d[:δo] = d[:δ0]

	d[:αo] = d[:γo] + d[:δo]
	d[:αu] = d[:γu] + d[:δu]

	d[:β]   = d[:βuo] = d[:βou] = d[:βoo] = d[:βuu]     # under 65 / over 65 infection rate
	#d[:αu]   = d[:γ0] * d[:ρ0] # Recovery rate + Death rate for old people
	#ay=gy+Dy  #Recovery rate + Death rate for young people
	d[:α0]   = d[:αu]*d[:ρu] + d[:αo]*d[:ρo]
	d[:go] = d[:gu] = 1. / d[:recovery_time] #recovery rate
    d[:No0] = d[:N] * d[:ρo]  # initial population over 65   (called No0 and Ny0)
    d[:Nu0] = d[:N] * d[:ρu]  # initial population under 65

	d[:βoo] = d[:mult_init] * (d[:αo]*d[:ρo] +  d[:αu]*d[:ρu]) / d[:Sf]
	d[:βuu] = d[:mult_init] * (d[:αu]*d[:ρu] +  d[:αu]*d[:ρu]) / d[:Sf]
	d[:βou] = d[:βuo] = min(d[:βuu], d[:βoo])

	#Initial Conditions. But I will call the initial condition function later,
	#which will use these paramters
	d[:Mo0] = 0                  # number of dead at t=0
	d[:Mu0] = 0
	d[:So0] = d[:No0] * d[:Sf]   # Total number of susceptibles at t=0
	d[:Su0] = d[:Nu0] * d[:Sf]
	d[:Ro0] = d[:No0] - d[:So0]  # Total number of recovered at =t0
	d[:Ru0] = d[:Nu0] - d[:Su0]
	#Initial fraction recovered population
	#d[:Io0] = d[:No0] - d[:So0] - d[:Ro0] - d[:Mo0]
	#d[:Iu0] = d[:Nu0] - d[:Su0] - d[:Ru0] - d[:Mu0]
	d[:Io0] = 1000
	d[:Iu0] = 1000
	# Make sure that the number of infected + number susceptible = total population at t=0, d[:N]
	d[:So0] -= d[:Io0]
	d[:Su0] -= d[:Iu0]
	#println("globalParameters, sum (I+S+R+M) old+young: ",
	#d[:Io0]+d[:Iu0]+d[:Ro0]+d[:Ru0]+d[:So0]+d[:Su0]+d[:Mo0]+d[:Mu0])
end

function initialConditions(p::Dict{Symbol,Float64})
	 # args p are parameters
     I0 = p[:Io0]
     S0 = p[:ρo] * p[:N]
     R0 = p[:Ro0]
     M0 = p[:Mo0]

	 i0 = p[:Iu0]
     s0 = p[:ρu] * p[:N]
     r0 = p[:Ru0]
     m0 = p[:Mu0]
     uinit = vcat(I0,S0,R0,M0,i0,s0,r0,m0)
     return uinit
end


function Covid!(du, u, p, t)
    #@timeit to "Covid!" begin # Must be at least one space or won't work
	I,S,R,M,i,s,r,m = u
    du[1] =  (p[:βuo]*S*i + p[:βoo]*S*I)/p[:N] - p[:αo]*I   # I
    du[2] = -(p[:βuo]*S*i + p[:βoo]*S*I)/p[:N]              # S
    du[3] = +p[:γo]*I                                           # R
    du[4] = +p[:δo]*I                                           # M
    du[5] =  (p[:βuu]*s*i + p[:βuo]*s*I)/p[:N] - p[:αu]*I   # i
    du[6] = -(p[:βuu]*s*i + p[:βuo]*s*I)/p[:N]              # s
    du[7] = +p[:γu]*I                                           # r
    du[8] = +p[:δu]*I                                           # m
    #du[8] = p[:N] - sum(du[2:8])
    #old   = I+S+R+M
    #young = i+s+r+m

      #d/dt [i + s + r + m] = (-p[:αu]*I + p[:γu] + p[:δu]) * I = 0
      #d/dt [I + S + R + M] = (-p[:αo]*I + p[:γo] + p[:δo]) * I = 0
      #.
      #Old people stay old; young people stay young, so the balance
      #.   above holds for each class individually.
      #Therefore, αo = γo + δo
      #.   and     αu = γu + δu
      #.
      #.    I + S + R + M = const = p[:ρo]
      #.    i + s + r + m = const = p[:ρu]
end

pp = globalParameters()
tspan = (0., 355.)

function solveProblem(tspan, u0, pp)
	prob = ODEProblem(Covid!, u0, tspan, pp)  # d is a parameter dictionary
	@time sol = solve(prob, Tsit5())
	return sol
end

for boo in .1:.01:.8
    pp = globalParameters()
    pp[:βuu] = boo
    # the cost of reinitializing the parameters is much lower than the solution cost
    derivedParameters(pp)
    u0 = initialConditions(pp)
    solveProblem(tspan, u0, pp)
	println("boo= ", boo)
end

function plotSolution(sol)
	labels=[:I,:S,:R,:M,:i,:s,:r,:m]
	p = plot()
	for (i,label) ∈ enumerate(labels)
		#print(i)
		plot!(p, sol.t, sol[i,:], label=label)
	end
	p
end

sol = solveProblem(tspan, u0, pp)
p = plotSolution(sol)
