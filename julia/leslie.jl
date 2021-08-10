using LinearAlgebra
using Plots
using RData


function make_leslies(sx_m, sx_f, asfr, srb, f_idx)
  # Average fertility between years, applying survival
  asfr_bar = (asfr .+ vcat(sx_f[f_idx[1:end - 1], :] .* asfr[1:end - 1, :], repeat([0.0], size(asfr)[2])')) .* 0.5
  f_bar_m = (1.0 .- (1.0 .+ srb).^(-1))' .* asfr_bar
  f_bar_f = ((1.0 .+ srb).^(-1))' .* asfr_bar
  i = 1
  make_leslie(f_bar_m[:,i], f_bar_f[:,i], sx_m[:,i], sx_f[:,i])
  leslies = [make_leslie(f_bar_m[:,i], f_bar_f[:,i], sx_m[:,i], sx_f[:,i]) for i in 1:size(sx_m)[2]]

  return leslies
end

function make_leslie(f_bar_m, f_bar_f, sx_m, sx_f)
  age_n = size(sx_m)[1]

  # Construct blocks of the Leslie matrix
  A = Matrix(vcat(
      transpose(repeat([0.0], size(sx_m)[1])),
      hcat(
          Diagonal(sx_m[1:end - 1]),
          vcat(repeat([0.0], size(sx_m)[1] - 2), sx_m[end])
      )
  ))
  A[1, f_idx] = f_bar_m
  B = zeros(age_n, age_n)
  C = copy(B)
  C[1,f_idx] .= f_bar_f
  D = Matrix(vcat(
      transpose(repeat([0.0], age_n)),
      hcat(
          Diagonal(sx_f[1:end - 1]),
          vcat(repeat([0.0], age_n - 2), sx_f[end])
      )
  ))
  leslie = vcat(
      hcat(A, B),
      hcat(C, D)
  )
  return leslie
end


data = load("data.rda", convert = true)
mod_data = data["temp"]

sx_f = mod_data["Sx"][:,1, :]
sx_m = mod_data["Sx"][:,2, :]
asfr = mod_data["asfr"]
srb = mod_data["srb"][1, :] ./ mod_data["srb"][2, :]
f_idx = mod_data["ss"]["p.fert.idx"]

leslies = make_leslies(sx_m, sx_f, asfr, srb, f_idx)

x =  vcat(mod_data["basepop"][:,1], mod_data["basepop"][:,2])

function project_pop(x, leslies)
  N = size(leslies)[1]
  X = Array{Float64}(undef, length(x), N)
  X[:, 1] .= x
  for i = 1:N
    x .= L * x
    X[:, i] .= x
  end
  return X
end
@time X = project_pop(x, leslies)

x_range = collect(mod_data["ss"]["proj_start"]:(mod_data["ss"]["proj_start"] + mod_data["ss"]["PROJ_YEARS"] - 1))
y_range = vec([string(i, j) for i = collect(15:80), j = ["m", "f"]])
Plots.heatmap(x_range, y_range, log10.(X))
heatmap(X[1:Int(size(X)[1] / 2), :] ./ X[Int(size(X)[1] / 2 + 1):end, :])

using Interpolations
A_x1 = 1:1.:10
A_x2 = 1:.5:20
g(x1, x2) = log(x1+x2)
A = [g(x1,x2) for x1 in A_x1, x2 in A_x2]
itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
sitp = scale(itp, A_x1, A_x2)
sitp(5., 10.) # exactly log(5 + 10)
sitp(5.6, 7.1) # approximately log(5.6 + 7.1)
sitp(5.6111, 7.0234)
log(5.6111 + 7.0234)
heatmap(A)
surface(A)
