using DifferentialEquations, Plots, BenchmarkTools
include("GlobalParameters.jl")
import .globalParameters_mod: globalParametersCorrect, globalParametersTest
getParams = globalParametersTest
#---------------------------------------------------------

function CIh!(du, u, p, t)
    #C, Ce, I, h = u
    C  = @view(u[1:1:2])
    Ce = @view(u[3:1:4])
    I  = @view(u[5:1:6])
    h  = @view(u[7:1:8])
    du[1:2] = [0,0]
    du[3:4] = [0,0]
    du[5:6] = [0,0]
    du[7:8] = [0,0]
end

nb_neurons = 2
eqs_per_neuron = 4
# Initial Conditions
c0 = rand(nb_neurons)
ce0 = rand(nb_neurons)
h0 = rand(nb_neurons)
I0 = rand(nb_neurons)
# Initial condition vector
u0 = vcat(c0,ce0,h0,I0)

pars = [0]
tspan = (0., 1.)  # make these real
prob = ODEProblem(CIh!, u0, tspan, pars)  # d is a parameter dictionary
println(prob)
@time sol = solve(prob, Tsit5())
println(sol.u[1])
println(sol.u[2])
println(sol.u[3])
for i ∈ 1:8
    println(i, sol[i,:])
end
println("--------------------")
plot(sol[:,1]) # u[1](t)
