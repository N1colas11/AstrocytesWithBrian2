using DifferejntialEquations, Plots, BenchmarkTools
include("GlobalParameters.jl")
import .globalParameters_mod: globalParametersCorrect, globalParametersTest
getParams = globalParametersTest
#---------------------------------------------------------
function coupled!(du, u, p, t)
    α, β, D1, D2 = p
    c1, c2 = u
    dc1  = α*c1 + D1 * (c2-c1)
    dc2 = β*c2  + D2 * (c1-c2)
    du .= [dc1, dc2]  # note the broadcast operator!
end

var = []
function coupledd!(du, u, p, t)
    α, β, D1, D2 = p
    c1, c2 = u
    dc1  = α*c1 + D1 * (c2-c1)
    dc2 = β*c2  + D2 * (c1-c2)
    push!(var, c2)
    du .= [dc1, dc2]  # note the broadcast operator!
end

u0 = rand(2)
tspan = (0.,10.)
pars = [-.2, -.1, 1., 1.]
prob = ODEProblem(coupledd!, u0, tspan, pars)
sol = solve(prob, Tsit5())
println("size(l_dc2)= ", size(l_dc2))
println("nb time steps= ", size(sol.t)








function Hill(x, p, n)
    return @. x^n / (x^n + p^n)
end

lJp = []
lJr = []
lJ1 = []
lCount = []

function CIh!(du, u, p, t)
    #C, Ce, I, h = u
    C = @view(u[1:1:2])
    Ce = @view(u[3:1:4])
    I = @view(u[5:1:6])
    h = @view(u[7:1:8])
    push!(lCount, 1)
    #println("... Ce= ", Ce)

    Tot_C = C
    Tot_CE = Ce               #: mole/meter**3  (sum of Ce in two half compartments in Brian2 code)

    # C: calcium
    # Ce: calcium in ER
    # h : fraction of open channels
    # h : IP3 concentration

    # rhs_C, etc. Are arrays, created via implicit vectorization involving numpy arrays
    dR2 = p[:R]^2 - p[:r]^2  #: meter**2
    s = p[:R] + p[:r]
    V = @. p[:Vrest] + (p[:F]*s/p[:Cm]) * (C - p[:Crest]) #: volt
    VVT = -2. * V/p[:V_T] #: 1  # value of 0, which leads to a zero denominator

    # Params:  r, P_r, p_open, dv_ER, P_CE, r, dv_ER, N_A, Omega_u, eta_p, K_P
    Jr   = @. (2*p[:r]/dR2) * p[:P_r] * p[:p_open] * (Ce - C) #: mole/meter**3/second  # p_open, Pr
    xxx  = @. -2. * p[:dv_ER]/VVT
    J1   = @. (4 *p[:P_CE]/p[:r]) * VVT * p[:dv_ER] * (C * exp(-2 * p[:dv_ER] / VVT) - Ce) / (exp(-xxx) - 1.) #: mole/meter**3/second  # dv_ER
    Jp   = @. (2 *p[:r]*p[:d_ER])/(p[:N_A]*dR2) * p[:Omega_u] * p[:eta_p] * C^2 / (C^2 + p[:K_P]^2) #: mole/meter**3/second # d_ER, N_A, Omega_u, eta_p, K_p

    # All three currents have the same units
    # (r/R^2) P_r p_open * C = (P_CE/r) * Volt * C  = (r * ER) (N_A * r^2) * Omega_u * eta_p

    push!(lJr, Jr)
    push!(lJp, Jp)
    push!(lJ1, J1)

    #println("Jr= ", Jr)
    #println("J1= ", J1)
    #println("Jp= ", Jp)
    rhs_C = @. Jr + J1 - Jp  #: mole / meter**3

    # Endoplasmic Reticulum Dynamnics
    # J1 is electro coupling (turn off because of exponential)
    rhs_Ce = @. Jp/(p[:rho]*p[:Lambda]) - (Jr + 0. * J1)  #:
    Jbeta  = 0  #*mole/meter**3/second  : mole/meter**3/second
    Jdelta = @. p[:o_delta] * (p[:k_delta]/(I+p[:k_delta])) * (C^2/(C^2+p[:K_delta]^2)) #: mole/second  # not sure about units, o_delta, k_delta, K_delta
    J5P    = @. p[:o_5P] * (I/(I+p[:K_5P])) #: mole/second # o_5P, K_5P
    J3K    = @. p[:o_3K] * (C^4/(C^4+p[:K_D]^4)) * (I/(I+p[:K_3])) #: mole/second # o_3K, K_D, K_3

    # assume L constant (CHECK)
    Tot_I  = I                #: mole/meter**3  (sum of Ce in two half compartments in Brian2 code)
    coupling_I  = @. (4*p[:D_I]  / p[:L]^2) * (I  - I)  # : mole/second/meter**3 (summed) (CHECK)

    # Open Channel dynamics
    OmegaH = @. (p[:Omega_2]*(I+p[:d_1]) + p[:O_2]*(I+p[:d_3])*C) / (I + p[:d_3]) #: Hz # Omega_2, O_2, d_1, d_3
    hinf   = @. p[:d_2] * (I + p[:d_1]) / (p[:d_2]*(I + p[:d_1]) + (I + p[:d_3])*C) #: 1 # d_2, d_1, d_3

    #-----------------------------
    # Handle diffusion terms, which involve terms from two adjacent compartments
    A1 = 0
    B1 = 0
    C1 = 0

    for i ∈ [1,2]
        A1 += @. (0.5 * p[:F] * s)/(p[:Cm] * p[:V_T]) * (dR2 / p[:L])         #: meter**4 / mole (summed)
        # Default: fac=1. I think it should be 1000 or 0.001 because V_T is givin in mVolts
        fac = 1
        #println("C[i]= ", C[i])
        B1 += @. (1. - (fac * s * p[:F] * C[i])/(2. * p[:Cm] * p[:V_T]))            #: meter (summed)
        C1 += @. dR2 * C[i] / p[:L]                              #: mole / meter**2 (summed)
    end

    #print("A1,B1,C1= ", A1, B1, C1)

    nb_connections = 1
    CC0               = @. ((-B1 + sqrt(B1^2 + 4*C1*A1)) / (2*A1)) / nb_connections #: mole / meter**3 (summed)
    V0                = @. (p[:Vrest] + (CC0 - p[:Crest]) * (p[:F] * s) / p[:Cm]) / nb_connections   #: volt  (summed)
    coupling_Ce       = @. (4*p[:D_CE] / p[:L]^2) * (Ce[1]-Ce[2]) * (1,-1)   #: mole/second/meter**3 (summed)  (CHECK) (most be Tot_Ce from each compartment
    electro_diffusion = @. -p[:P_Ca] * V / (p[:R] * p[:V_T]) * (Ce*exp(-2*V/p[:V_T]) - C) / (exp(-2*V/p[:V_T]) - 1) #: mole / second / meter**3
    coupling_electro  = @. (4*p[:D_C]/p[:L]^2/p[:V_T]) * (CC0 + C) * (V0 - V)  #: mole/second/meter**3 (summed)
    coupling_C        = @. (4*p[:D_C] / p[:L]^2) * (C[1] - C[2]) * (1,-1) #: mole/second/meter**3 (summed)  # ERROR. Cannot be zero.

    dh  = @. OmegaH * (hinf - h) #: 1
    dI  = @. (Jbeta + Jdelta - J3K - J5P) / (p[:Lambda]*(1-p[:rho])) #: mole/meter**3
    dC  = @. 0*coupling_C + 0*coupling_electro + 0*electro_diffusion + Jr + J1 - Jp  #: mole / meter**3
    dCE = @. 0*coupling_Ce + Jp/(p[:rho]*p[:Lambda]) - (Jr + J1)  #: mole/meter**3
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
    dC  = @. 1. *Jr + 1. *J1 - 1. *Jp / p[:Lambda]
    dCE = @. 0. *coupling_Ce + 1. *Jp/(p[:rho]*p[:Lambda]) - 1. * (Jr + J1)  #: mole/meter**3
    #  Jbeta: always zero
    #  Jdelta: goes to 3^4
    #  J3K: ``goes to -10^-12  (NOT GOOD)
    #  J5P: ``goes to -10^-13  (NOT GOOD) (goes to -.6 if I do not normalized by volume)

    println("Jbeta, Jdelta, J3K, J5P = ", Jbeta, Jdelta, J3K, J5P)
    dI  = @. (0. *Jbeta + 0. *Jdelta - 0. *J3K - 1. *J5P) # / (p[:Lambda]*(1-p[:rho])) #: mole/meter**3
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

u0 = initialConditions()
pars = getParams()
tspan = (0., 2.e-7)  # make these real
#u0 = [.1, .2, .3, .4]``
# We have 8 equations. How to collect them?
#CIh!(du, u0, p, t)
prob = ODEProblem(CIh!, u0, tspan, pars)  # d is a parameter dictionary
#@time sol = solve(prob, Tsit5())
@time sol = solve(prob, Euler(), dt=1.e-7)
println(prob)
#interp = tspan[1]:(tspan[2]-tspan[1])/10.:tspan[2]
#sol_interp = sol(interp)
println("t= ", sol.t)
println("--------------------")
# Not using sol_interp correctly!
#plot(sol_interp[1,:])
p1 = plot(sol.t, sol[1,:], label="C")
p2 = plot(sol.t, sol[3,:], label="Ce")
p3 = plot(sol.t, sol[5,:], label="I")
p4 = plot(sol.t, sol[7,:], label="h")
plot(p1,p2,p3,p4, layout=(2,2))
#plot!(p, sol.t, sol[7,:])
print(sol.t)
print(lJ1[1,:])
print(size(lJ1[:,1]))
plot(sol.t, lJ1)

#for i ∈ 1:8
    #print("i= ", i, ", sol= ", sol[i,:], "\n")
#end
#
function solve_CIh!(u0, tspan, d)
    solve(prob, Tsit5())
end




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

run_problem(αs, βs)

t = [-1,0,1]
a = test_plot()
plot(t, a[1])
plot!(t, a[2])

plot(a[2])

ENV["GRDIR"]=""
Pkg.build("GR")

a = [1,2,3]
b = [2,3,4]
println(@. 3 + b + 5 * log(a))
