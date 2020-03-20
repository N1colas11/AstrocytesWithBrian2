module globalFunctions_mod

using Plots
using Random
using BenchmarkTools
using Profile, ProfileView
using Distributions

function plotSol(sol)
	p1 = plot(sol.t, sol[1,:], label="C1")
	#p1 = plot!(p1, sol.t, sol[2,:], label="C2")
	p2 = plot(sol.t, sol[3,:], label="Ce1")
	#p2 = plot!(sol.t, sol[4,:], label="Ce2")
	p3 = plot(sol.t, sol[5,:], label="I1")
	p3 = plot!(p3,sol.t, sol[6,:], label="I2")
	p4 = plot(sol.t, sol[7,:], label="h1")
	p4 = plot!(p4,sol.t, sol[8,:], label="h2")
	plot(p1,p2,p3,p4, layout=(2,2))
end
#---------------------------------------
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
#----------------------------------------------------------------------
#=
function PoissonPulses1(t::Float64, half_width::Float64)
    REAL t_s
    Branch& b = *current_branch # ?????

    # looping through values and returning 1. if t is half_width within a Poisson event time
    # much faster than summing Poisson Pulses
    sum::Float64 = 0.
    for i in 1:length(poisson_times1)
        #if 0
            if (fabs(b.poisson_times1[i] - t) < half_width)
                #("half_width=%f\n",half_width);
                #("%s t_k=%f and t=%f, so the difference: %f < %f\n", b.name.c_str(), b.poisson_times1[i], t, fabs(b.poisson_times1[i] - t), half_width );
                return 1.
			end
        # else
            sum = sum + ( 1./(fabs(half_width) * sqrt(π)) ) * exp( -pow(t/half_width,2.) )
        #endif
   end
   return 0.0 # if you make it through the for loop there are no events
end
=#
#---------------------------------------------------------------
function generatePoissonEvents!(tk::Vector{Float64}, λ::Float64)
	# Generate time spikes according to a Poisson Process
	# 1) generate nb_events in U(0.,1.)
	# For each of nb_times time intervals, check the next random number against λ δt and create an
	#   event if λ δt < random number

	# Input: tk: a preallocated vector of size equal to the number of desired events
	# Later, I should generate them on demand using Iterators perhaps

	# Version 1. We will run at constant time steps.
	#println("rands: ", rands)

	# Get samples from an exponential distribution
	println("********* Gordon ***************")
	nb_events = length(tk)
	if nb_events == 0
		println("nb_events should be > 0")
		return -1
	end
	dist = Distributions.Exponential(λ)
	expon_samples = rand(dist, nb_events)
	shuffled_intervals = shuffle(expon_samples)

	# Initial time is zer0. List of Float64 of spike times.
	cumsum!(tk, shuffled_intervals)
	#println("Y= ", shuffled_intervals[1:10])
	#println("tk= ", tk[1:10])

	#plot tk (Must be a better way)
	pl = scatter(tk, x->1, title="Poisson train tk")
	display(pl)

	return 
end

#tk = Vector{Float64}(undef, 50)
#generatePoissonEvents!(tk, 3.)
#----------------------------------------------------------------------
function generatePoissonEvents(nb_events::Int64, λ::Float64, dt::Float64)
	# Generate time spikes according to a Poisson Process
	# 1) generate nb_events in U(0.,1.)
	# For each of nb_times time intervals, check the next random number against λ δt and create an
	#   event if λ δt < random number
	# dt: a fixed time step
	#dist = Distributions.Poisson(10.)
	#dist = Distributions.Normal(11., 3.)

	# Version 1. We will run at constant time steps.
	nb_trials = 20
	println("nb_trials: ", nb_trials)
	#trials = rand(dist, 10000)
	# Uniform distribution
	nb_rands = 100
	rands = rand(nb_rands)
	#println("rands: ", rands)
	R = λ * dt
	println("R= ", R)
	bools = @. rands < R
	#println("bools= ", bools)
	max_true = sum(bools)
	frac_true = max_true / nb_rands
	println("frac_true= ", frac_true)

	# get samples from an exponential distribution
	dist = Distributions.Exponential(λ)
	expon_samples = rand(dist, nb_rands)
	#println("expon_samples= ", expon_samples)

	# Randominze the order of the expon_samples (ISIs)
	#idx = randperm(expon_samples)
	shuffled_intervals = shuffle(expon_samples)
	#println("sum(Y1),Sum(Y2)= ", sum(Y1), ",  ", sum(Y2))

	# Compute Spike times, starting from zero
	bins = Vector{Float64}(undef,nb_trials)

	t0 = 0.
	tk = cumsum(shuffled_intervals)
	println("Y= ", shuffled_intervals[1:10])
	println("tk= ", tk[1:10])

	#plot tk (Must be a better way)
	y = [1. for i in 1:length(tk)]
	#pl = scatter(tk, y, title="Poisson train tk")
	pl = scatter(tk, x->1, title="Poisson train tk")
	display(pl)
	println(typeof(tk))

	return
end

#for i in 1:1
	#generatePoissonEvents(10, 2., .1)
#end
#---------------------------------------------------------------
function generatePoissonEvents!(
	poisson_times1::Vector{Float64}, poisson_times2::Vector{Float64},
		lambda1::Float64, lambda2::Float64, T::Float64, seed::Int64)
#=
******************************************************************************/
  Purpose:
    homo_poisson returns list of event times for n_inputs compartments given time interval T
    ***SEED RESETS EVERY FUNCTION CALL ***

  Discussion:
    The homogeneous poisson process used to generate event times
    (with intensity lamda)
  =#

#  RVEC poisson_times1;
  # generate poisson_times2
  #count::Int64 = 0
  t0::Float64  = 0.
  go           = true
  tn::Float64  = 0.

  # generate poisson_times1
  while true
	buffer = rand() #normal.r4_uniform_01 ( seed ); <<<<
	#("poisson buffer: %f",buffer);
	w      = -log(buffer) / lambda1

	if length(poisson_times1) > 1
		tnp1 = tn + w
	else
		tnp1 = t0 + w
	end

	# make sure new event time is < T
	if (tnp1 > T)
		break
 	else
		#count += 1
		#poisson_times1.push_back(tnp1)
		push!(poisson_times1, tnp1)
		#println("tnp1= ", tnp1)
		tn = tnp1 + 1
	end
  end

  #println("\n\npoisson_times1 = [")
  #for i in 1:length(poisson_times1)
 	#println("%f, ",poisson_times1[i])
  #end
  #println("]\n\n")


  while true
	buffer = rand()
	#("poisson buffer: %f",buffer)
	w = - log(buffer) / lambda2

	if length(poisson_times2) > 1
		tnp1 = tn + w
	else
		tnp1 = t0 + w
	end

	# make sure new event time is < T
	if tnp1 > T
		break
 	else
		#count += 1
		push!(poisson_times2, tnp1)
		tn = tnp1 + 1
	end
  end

  #println("\n\npoisson_times2 = [")
  #for (int i=0;i<poisson_times2.size();i++)
  #for i in 1:length(poisson_times2)
 	#println("%f, ",poisson_times2[i])
  #end
  #println("]\n\n")
  return 0.
end

#----------------------------------------------------------------------
function testPoisson()
	λ_1 = 2.
	λ_2 = .3
	#buffer1 = Vector{Float64}[]
	#buffer2 = Vector{Float64}[]
	buffer1 = Float64[]
	buffer2 = Float64[]
	T = 3.
	println(typeof(T))
	seed = 1000
	generatePoissonEvents!(buffer1, buffer2, λ_1 ,λ_2,  T, seed)
	#println("buffer1= ", buffer1)
	#println("buffer2= ", buffer2)
	return λ_1, λ_2, buffer1, buffer2
end

#@profiler for i in 1:10000
	#λ_1, λ_2, buffer1, buffer2 = testPoisson()
	##println("-------i= ", i, " ------")
#end


# end module
end
