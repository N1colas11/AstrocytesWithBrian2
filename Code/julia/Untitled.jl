using DifferentialEquations


function diffusion!(du, u, p, t)
    α = p
    du = -α*u
end

u0 = 3.0
tspan = (0., 2.)
p = 1.0

prob = ODEProblem(diffusion!, u0, tspan, p)
sol = solve(prob, Tsit5())
