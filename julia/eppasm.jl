module eppasm
export simmodJ, pulldata

using LinearAlgebra
using Serialization
using StaticArrays
using YAML

include("pop_project.jl")
include("hiv_mod.jl")
include("treatment.jl")

function simmodJ(data = nothing)::Dict
  if data === nothing
    data = pulldata()
  end
  p = prepp(data)
  out_dict = runsim!(p)
  return out_dict
end

function pulldata()
  data = deserialize("data")
  return data
end

function prepp(data)
  p = NamedTuple()
  # Global parameters
  gpar = YAML.load_file("julia/gpar.yml"; dicttype=Dict{Symbol,Any})
  for k = keys(gpar)
    p = merge(p, [k=>gpar[k]])
  end
  # Passed parameters
  for k = keys(data)
    if isa(data[k], Array)
      if ndims(data[k]) > 1
        p = merge(p, @SArray [k=>permutedims(data[k], reverse(1:ndims(data[k])))])
      else
        p = merge(p, @SArray [k=>data[k]])
      end
    else
      p = merge(p, [k=>data[k]])
    end
  end
  p = merge(p, [:PROJ_YEARS => p.ss[:PROJ_YEARS]])
  p = merge(p, [:HIVSTEPS_PER_YEAR => p.ss[:hiv_steps_per_year]])
  p = merge(p, [:hAG_SPAN => p.ss[:h_ag_span]])
  hAG_START = zeros(Int, p.hAG)
  for ha = 2:p.hAG
    hAG_START[ha] = hAG_START[ha-1] + p.hAG_SPAN[ha-1]
  end
  p = merge(p, [:hAG_START => hAG_START])
  p = merge(p, [:DT => 1.0 / p.HIVSTEPS_PER_YEAR])
  p = merge(p, [:everARTelig_idx => p.hDS])
  
  # Construct age binners
  age_binner = zeros(sum(p.hAG_SPAN), length(p.hAG_SPAN))
  max_age_binner = zeros(sum(p.hAG_SPAN), length(p.hAG_SPAN))
  j_start = 1
  for i = 1:length(p.hAG_SPAN)
    j_n = p.hAG_SPAN[i]
    for j = j_start:(j_n + j_start - 1)
      age_binner[j, i] = 1
      if j == (j_n + j_start - 1)
        max_age_binner[j, i] = 1
      end
    end
    j_start += j_n
  end
  p = merge(p, [:age_binner => age_binner])
  p = merge(p, [:max_age_binner => max_age_binner])
  p = merge(p, [:idx_15_49 => p.pIDX_15TO49:(p.pIDX_15TO49 + p.pAG_15TO49 - 1)])
  p = merge(p, [:h_idx_15_49 => (p.hIDX_15TO49 + 1):(p.hIDX_15TO49 + p.hAG_15TO49 + 1)])
  p = merge(p, [:fert_idx => p.pIDX_FERT:(p.pIDX_FERT + p.pAG_FERT - 1)])
  p = merge(p, [:h_fert_idx => (p.hIDX_FERT + 1):(p.hIDX_FERT + p.hAG_FERT)])
  return p
end

function runsim!(p)
  # Setup baseline
  x, out_dict = setupbaseline(p)
  
  # Simulate all years
  simall!(p, x, out_dict)

  return out_dict
end

function setupbaseline(p)
  # States
  pop = @MArray zeros(p.pDS, p.NG, p.pAG)
  hivpop = @MArray zeros(p.NG, p.hAG, p.hDS)
  artpop = @MArray zeros(p.NG, p.hAG, p.hDS, p.hTS)
  grad = @MArray zeros(p.NG, p.hAG, p.hDS)
  gradART = @MArray zeros(p.NG, p.hAG, p.hDS, p.hTS)

  # Fill base pop
  pop[p.HIVN, :, :] = p.basepop
  pop[p.HIVP, :, :] .= 0.

  x = Dict{Symbol, MArray}(
    :pop => pop,
    :hivpop => hivpop,
    :artpop => artpop,
    :grad => grad,
    :gradART => gradART
  )  

  # Outputs
  prev15to49::Array{Float64, 1} = zeros(p.PROJ_YEARS)
  pop_ts::Array{Float64, 4} = zeros(p.PROJ_YEARS, p.pDS, p.NG, p.pAG)
  hivpop_ts::Array{Float64, 4} = zeros(p.PROJ_YEARS, p.NG, p.hAG, p.hDS)
  artpop_ts::Array{Float64, 5} = zeros(p.PROJ_YEARS, p.NG, p.hAG, p.hDS, p.hTS)
  infections::Array{Float64, 3} = zeros(p.PROJ_YEARS, p.NG, p.pAG)
  prev15to49_ts::Array{Float64, 1} = zeros((p.PROJ_YEARS-1) * p.HIVSTEPS_PER_YEAR)
  incrate15to49_ts::Array{Float64, 1} = zeros((p.PROJ_YEARS-1) * p.HIVSTEPS_PER_YEAR)
  hivdeaths::Array{Float64, 3} = zeros(p.PROJ_YEARS, p.NG, p.pAG)
  natdeaths::Array{Float64, 3} = zeros(p.PROJ_YEARS, p.NG, p.pAG)
  incid15to49::Array{Float64, 1} = zeros(p.PROJ_YEARS)
  aidsdeaths_noart::Array{Float64, 4} = zeros(p.PROJ_YEARS, p.NG, p.hAG, p.hDS)
  aidsdeaths_art::Array{Float64, 5} = zeros(p.PROJ_YEARS, p.NG, p.hAG, p.hDS, p.hTS)
  artinit::Array{Float64, 4} = zeros(p.PROJ_YEARS, p.NG, p.hAG, p.hDS)
  pregprevlag::Array{Float64, 1} = zeros(p.PROJ_YEARS)
  pregprev::Array{Float64, 1} = zeros(p.PROJ_YEARS)
  popadjust::Array{Float64, 3} = zeros(p.PROJ_YEARS, p.NG, p.pAG)
  entrantprev_out::Array{Float64, 1} = zeros(p.PROJ_YEARS)
  
  # Fill with initial population
  pop_ts[1, :, :, :] = pop
  
  out_dict = Dict{Symbol, Array}(
    :hivpop => hivpop_ts,
    :artpop => artpop_ts,
    :infections => infections,
    :hivdeaths => hivdeaths,
    :natdeaths => natdeaths,
    :aidsdeaths_noart => aidsdeaths_noart,
    :aidsdeaths_art => aidsdeaths_art,
    :popadjust => popadjust,
    :artinit => artinit,
    :pregprevlag => pregprevlag,
    :incrate15to49_ts => incrate15to49_ts,
    :prev15to49_ts => prev15to49_ts,
    :rvec => p.rvec,
    :prev15to49 => prev15to49,
    :pregprev => pregprev,
    :incid15to49 => incid15to49,
    :entrantprev_out => entrantprev_out,
    :pop => pop_ts
  )

  return x, out_dict
end

function simall!(p, x, out_dict)
  for t = 2:p.SIM_YEARS
    simone!(t, p, x, out_dict)
    @assert sum(x[:pop][p.HIVP, :, :]) - (sum(x[:artpop]) + sum(x[:hivpop])) < 1
  end

  # Calculate additional outputs
  finaloutputs(p, out_dict)
end

function simone!(t, p, x, out_dict)
  # Project forward the population
  births_by_ha = pop_project_one_step!(p, t, x, out_dict)
  
  # HIV model simulation
  if p.proj_steps[(t - 2) * p.HIVSTEPS_PER_YEAR + 1] >= p.tsEpidemicStart
    hiv_mod_all(p, x, t, births_by_ha, out_dict)
  end
  
  # Adjust to a target population
  if p.popadjust
    adjust_pop!(p, t, x, out_dict)
  end

  # Prevalence in pregnant women
  preg_women_prev_one_step!(p, t, x, births_by_ha, out_dict)

  # Update output timeseries objects
  out_dict[:pop][t, :, :, :] = x[:pop]
  out_dict[:hivpop][t, :, :, :] = x[:hivpop]
  out_dict[:artpop][t, :, :, :, :] = x[:artpop]
end

function finaloutputs(p, out_dict)
  hivn15to49 = dropdims(sum(out_dict[:pop][:, p.HIVN, :, p.pIDX_15TO49:(p.pIDX_15TO49 + p.pAG_15TO49 - 1)], dims = (2, 3, 4)), dims = (2, 3))
  hivp15to49 = dropdims(sum(out_dict[:pop][:, p.HIVP, :, p.pIDX_15TO49:(p.pIDX_15TO49 + p.pAG_15TO49 - 1)], dims = (2, 3, 4)), dims = (2, 3))
  out_dict[:prev15to49] = hivp15to49 ./ (hivn15to49 + hivp15to49)
  out_dict[:incid15to49][2:end] ./= hivn15to49[1:end - 1]
end

end
