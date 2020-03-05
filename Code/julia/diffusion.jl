using DifferentialEquations
using Plots
using BenchmarkTools

function diffusion!(du, u, p, t)
    α, β = p
    c, ce = u
    dc  = -α*c
    dce = -β*ce
    du .= [dc, dce]  # note the broadcast operator!
end

function coupled!(du, u, p, t)
    α, β, D1, D2 = p
    c1, c2 = u
    dc1  = α*c1 + D1 * (c2-c1)
    dc2 = β*c2  + D2 * (c1-c2)
    du .= [dc1, dc2]  # note the broadcast operator!
end

u0 = [1., 2.]
tspan = (0., 3.)
p = [1.0, 1.]
p1 = [1., 2., 8., 5.]  # α, β, D1, D2

prob = ODEProblem(coupled!, u0, tspan, p1)
sol = solve(prob, Tsit5())
plot(sol)
