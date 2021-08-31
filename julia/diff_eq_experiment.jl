using DifferentialEquations
using DiffEqSensitivity
using LinearAlgebra
using Serialization
using Plots
using ModelingToolkit

data = deserialize("data")

function f!(du, u, p, t)
  # Aging, dying, and migration
  du[1, :] .+= p[:entrantpop][:, t + 1]
  du[1:end - 1, :] .-= u[1:end - 1, :]
  du[2:end, :] .+= u[1:end - 1, :] - u[1:end - 1, :] .* (1 .- p[:Sx][1:end - 1, :, t + 1])
  du[end, :] .-= u[end, :] .* (1 .- p[:Sx][end, :, t + 1])
  du .+= p[:netmigr][:, :, t + 1]
end
u0 = data[:basepop]
tspan = (0, 53)
prob = DiscreteProblem(f!, u0, tspan, data, jac = true)
@time sol = solve(prob)
heatmap(sol[:, 1, :])

@variables t u(t)
@parameters α β q
# D = Difference(t; dt = 1)
D = Differential(t)

@named popproj = ODESystem([
  D(u) ~ (-u + α)*q + β
])