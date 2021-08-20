module Tst
  include("eppasm.jl")
  using .eppasm
  using Plots
  using BenchmarkTools
  plot_diagnostics = true
  data = pulldata()
  @time simmodJ(data)
  
  data = pulldata()
  @time out = simmodJ(data)
  out[:infections][7, :, :]'
  out[:incid15to49]
  out[:rvec][25]
  out[:incrate15to49_ts][100]
  out[:popadjust][3, :, :]
  @time simmodJ()
  sum(out[:artpop], dims = (2, 3, 4, 5))[33]

  correct_prev = [0.00045, 0.00080, 0.0014, 0.00245, 0.00424, 0.00725, 0.01214,
  0.01985, 0.03147, 0.04804, 0.07013, 0.0975, 0.12857, 0.16083,
  0.19143, 0.21797, 0.23908, 0.25419, 0.26363, 0.26817, 0.26871,
  0.26618, 0.26129, 0.2551, 0.24887, 0.24365, 0.23928, 0.23545,
  0.23186, 0.22866, 0.22569, 0.22277, 0.21931, 0.21555, 0.21147,
  0.2069, 0.20155, 0.19559, 0.1893, 0.1829, 0.1763, 0.16953, 0.16266]

  # artpop = dropdims(sum(out[:artpop], dims = (2, 3,  4, 5)), dims = (2, 3, 4, 5))
  # println("Input:", round.(p[:art15plus_num], digits = 3))
  # println("Output:", round.(artpop, digits = 3))
# p = prep_p()
# p[:art15plus_num] .= 0
# @time out = simmodJ(p)

if plot_diagnostics
  plot(out[:prev15to49][11:53])
  vline!([data[:tARTstart]])
  display(plot!(correct_prev))
  plot_pop_ta = sum(out[:pop], dims = (2, 3))[:, 1, 1, :]
  display(heatmap(plot_pop_ta'))
end
# println(round.(out[:prev15to49], digits = 3))
artpop = dropdims(sum(out[:artpop], dims = (2, 3, 4, 5)), dims = (2, 3, 4, 5))
# hivpop = dropdims(sum(out[:hivpop], dims = (2, 4)), dims = (2, 4))
# heatmap(hivpop')
plot(artpop)
plot!(dropdims(sum(data[:art15plus_num] .* .!data[:art15plus_isperc], dims = 1), dims = 1))
display(plot!(dropdims(sum(data[:art15plus_num] .* data[:art15plus_isperc] .* dropdims(sum(out[:pop][:, 2, :, :], dims = 3), dims = 3)', dims = 1), dims = 1)))
end