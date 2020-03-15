using DifferentialEquations, ModelingToolkit

@parameters t σ ρ β
@variables x(t), y(t), z(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z) - y,
       D(z) ~ x*y - β*z]

eqs
