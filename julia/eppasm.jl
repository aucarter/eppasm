module eppasm
export simmodJ, prep_par

using LinearAlgebra
using RCall
using Serialization

include("pop_project.jl")
include("hiv_mod.jl")
include("treatment.jl")

const AGE_START = 15 # FIRST
const NG = 2 # Number of Genders
const pAG = 66 # Number of Ages
const pDS = 2 # Number of population HIV states (pos / neg)
const pIDX_FERT = 1 # Start fertility idx
const pAG_FERT = 35 # Number of fertile ages
const pIDX_15TO49 = 1 # Start 15:49 idx
const pAG_15TO49 = 35 # Number of 15:49 ages
const pIDX_15PLUS = 1 # Start 15+
const pAG_15PLUS = 66 # Number of 15+ AGES
const hAG = 9 # HIV sim age bins
const hDS = 7 # CD4 count bins
const hTS = 3 # Time since ART initiation
const hIDX_FERT = 0
const hAG_FERT = 8
const hIDX_15TO49 = 0
const hAG_15TO49 = 8
const hIDX_15PLUS = 1
const hAG_15PLUS = 9
const hIDX_CD4_350 = 3
const MALE = 1
const FEMALE = 2
const HIVN = 1
const HIVP = 2
const ART0MOS = 1
const EPP_RSPLINE = 0
const EPP_RTREND = 1
const EPP_DIRECTINCID = 2 # annual direct incidence inputs (as Spectrum)
const INCIDMOD_EPPSPEC = 0
const INCIDMOD_TRANSM = 1
const INCIDPOP_15TO49 = 0 # age range corresponding to incidence input
const INCIDPOP_15PLUS = 1
const CD4_LOW_LIM = [500, 350, 250, 200, 100, 50, 0]
const CD4_UPP_LIM = [1000, 500, 350, 250, 200, 100, 50]
const eppmod = 0
const incidmod = 0
const idx_15_49 = pIDX_15TO49:(pIDX_15TO49 + pAG_15TO49 - 1)
const h_idx_15_49 = (hIDX_15TO49 + 1):(hIDX_15TO49 + hAG_15TO49 + 1)

function prep_par()
  # par = rcopy(R"""
  # library(eppasm)
  # pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
  # bw <- prepare_spec_fit(pjnz, proj.end=2022.5)
  # bw_fp <- attr(bw$Urban, "specfp")
  # bw_fp <- prepare_rhybrid(bw_fp)
  # bw_fp$ancsitedata <- TRUE
  # bw_fp$ancrt <- "census"
  # bw_fp$ancrtsite.beta <- 0
  # bw_fp$logitiota <- TRUE
  # bw_theta <- c(-0.407503322169364, -2.76794181367538, -1.26018073624346, 1995.96447776502,
  #               -0.00307437171215574, 0.0114118307148102, 0.00760958379603691, 0.02,
  #               2.24103194827232, -0.0792123921862689, -5.01917961803606, 0.359444135205712,
  #               -6.10051517060137)
  # param <- fnCreateParam(bw_theta, bw_fp)
  # bw_fp <- update(bw_fp, list=param)
  # bw_fp
  # """)
  par  = deserialize("test")
  return par
end

function sim_one!(t, par, pop, hivpop, artpop,
  everARTelig_idx, age_binner, max_age_binner, grad, out_dict)

  # Project forward the population
  births, births_by_ha, last_hivpop, last_artpop = pop_project_one_step!(par, t, pop, hivpop, artpop,  age_binner, max_age_binner, out_dict)
  ## HIV model simulation
  cd4elig_idx = par[:artcd4elig_idx][t] 
  anyelig_idx = (par[:specpop_percelig][t] > 0 || par[:pw_artelig][t] > 0) ? 1 : (par[:who34percelig] > 0) ? hIDX_CD4_350 : cd4elig_idx
  everARTelig_idx = (anyelig_idx < everARTelig_idx) ? anyelig_idx : everARTelig_idx
  if par[:proj_steps][(t - 2) * HIVSTEPS_PER_YEAR + 1] >= par[:tsEpidemicStart]
    hiv_mod_all(par, t, anyelig_idx, cd4elig_idx, artpop, hivpop, everARTelig_idx, births_by_ha, pop, age_binner, grad, out_dict)
  end
  
  # Adjust to a target population
  if par[:popadjust]
    adjust_pop!(par, t, pop, hivpop, artpop, out_dict)
  end

  # Prevalence in pregnant women
  preg_women_prev_one_step!(par, t, pop, hivpop, artpop, births, last_hivpop, last_artpop, out_dict)

  # Update output timeseries objects
  out_dict[:pop][t, :, :, :] = pop
  out_dict[:hivpop][t, :, :, :] = hivpop
  out_dict[:artpop][t, :, :, :, :] = artpop
end

function setup_baseline(par)
  pop = zeros(pDS, NG, pAG)
  pop[HIVN, :, :] = par[:basepop]'
  pop[HIVP, :, :] .= 0.
  hivpop = zeros(NG, hAG, hDS)
  artpop = zeros(NG, hAG, hDS, hTS)
  prev15to49 = zeros(PROJ_YEARS)
  pop_ts::Array{Float64, 4} = zeros(PROJ_YEARS, pDS, NG, pAG)
  hivpop_ts::Array{Float64, 4} = zeros(PROJ_YEARS, NG, hAG, hDS)
  artpop_ts::Array{Float64, 5} = zeros(PROJ_YEARS, NG, hAG, hDS, hTS)
  infections::Array{Float64, 3} = zeros(PROJ_YEARS, NG, pAG)
  prev15to49_ts::Array{Float64, 1} = zeros((PROJ_YEARS-1) * HIVSTEPS_PER_YEAR)
  incrate15to49_ts::Array{Float64, 1} = zeros((PROJ_YEARS-1) * HIVSTEPS_PER_YEAR)
  hivdeaths::Array{Float64, 3} = zeros(PROJ_YEARS, NG, pAG)
  natdeaths::Array{Float64, 3} = zeros(PROJ_YEARS, NG, pAG)
  incid15to49::Array{Float64, 1} = zeros(PROJ_YEARS)
  aidsdeaths_noart::Array{Float64, 4} = zeros(PROJ_YEARS, NG, hAG, hDS)
  aidsdeaths_art::Array{Float64, 5} = zeros(PROJ_YEARS, NG, hAG, hDS, hTS)
  artinit::Array{Float64, 4} = zeros(PROJ_YEARS, NG, hAG, hDS)
  pregprevlag::Array{Float64, 1} = zeros(PROJ_YEARS)
  pregprev::Array{Float64, 1} = zeros(PROJ_YEARS)
  popadjust::Array{Float64, 3} = zeros(PROJ_YEARS, NG, pAG)
  entrantprev_out::Array{Float64, 1} = zeros(PROJ_YEARS)
  grad::Array{Float64, 3} = zeros(NG, hAG, hDS)
  everARTelig_idx = hDS

  # Fill with initial population
  pop_ts[1, :, :, :] = pop

  # Construct age binners
  age_binner = zeros(sum(hAG_SPAN), length(hAG_SPAN))
  max_age_binner = zeros(sum(hAG_SPAN), length(hAG_SPAN))
  j_start = 1
  for i = 1:length(hAG_SPAN)
    j_n = hAG_SPAN[i]
    for j = j_start:(j_n + j_start - 1)
      age_binner[j, i] = 1
      if j == (j_n + j_start - 1)
        max_age_binner[j, i] = 1
      end
    end
    j_start += j_n
  end

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
    :rvec => par[:rvec],
    :prev15to49 => prev15to49,
    :pregprev => pregprev,
    :incid15to49 => incid15to49,
    :entrantprev_out => entrantprev_out,
    :pop => pop_ts
  )
  return pop, hivpop, artpop, everARTelig_idx, age_binner, 
    max_age_binner, grad, out_dict
end

function sim_all!(par, pop, hivpop, artpop, everARTelig_idx, age_binner, max_age_binner, grad, out_dict)
  for t = 2:par[:SIM_YEARS]
    sim_one!(t, par, pop, hivpop, artpop,
      everARTelig_idx, age_binner, max_age_binner, grad, out_dict)
  end
end
function run_sim!(par)
  # Setup baseline
  pop, hivpop, artpop, everARTelig_idx, age_binner, max_age_binner, grad, out_dict= setup_baseline(par)
  # Simulate all years
  sim_all!(par, pop, hivpop, artpop, everARTelig_idx, age_binner, max_age_binner, grad, out_dict)

  # Calculate outputs
  hivn15to49 = dropdims(sum(out_dict[:pop][:, HIVN, :, pIDX_15TO49:(pIDX_15TO49 + pAG_15TO49 - 1)], dims = (2, 3, 4)), dims = (2, 3))
  hivp15to49 = dropdims(sum(out_dict[:pop][:, HIVP, :, pIDX_15TO49:(pIDX_15TO49 + pAG_15TO49 - 1)], dims = (2, 3, 4)), dims = (2, 3))
  out_dict[:prev15to49] = hivp15to49 ./ (hivn15to49 + hivp15to49)
  out_dict[:incid15to49][2:end] ./= hivn15to49[1:end - 1]

  return out_dict
end

function prep_inputs!(par)
  if(size(par[:Sx])[1] == 66)
    par[:Sx] = permutedims(par[:Sx], [3, 2, 1])
  end
  if(size(par[:netmigr])[1] == 66)
    par[:netmigr] = permutedims(par[:netmigr], [3, 2, 1])
  end
  if(size(par[:cd4_mort])[1] != 2)
    par[:cd4_mort] = permutedims(par[:cd4_mort], (3, 2, 1))
  end
  if(size(par[:cd4_prog])[1] != 2)
    par[:cd4_prog] = permutedims(par[:cd4_prog], (3, 2, 1))
  end
  if(size(par[:frr_cd4])[1] != 53)
    par[:frr_cd4] = permutedims(par[:frr_cd4], (3, 2, 1))
  end
  if(size(par[:frr_art])[1] != 53)
    par[:frr_art] = permutedims(par[:frr_art], (4, 3, 2, 1))
  end
  if(size(par[:incrr_age])[1] != 53)
    par[:incrr_age] = permutedims(par[:incrr_age], (3, 2, 1))
  end
  if(size(par[:cd4_initdist])[1] != 2)
    par[:cd4_initdist] = permutedims(par[:cd4_initdist], (3, 2, 1))
  end
  if(size(par[:targetpop])[1] == 66)
    par[:targetpop] = permutedims(par[:targetpop], [3, 2, 1])
  end
  if(size(par[:art_mort])[1] != 2)
    par[:art_mort] = permutedims(par[:art_mort], (4, 3, 2, 1))
  end
  if(size(par[:paedsurv_artcd4dist])[1] != 53)
    par[:paedsurv_artcd4dist] = permutedims(par[:paedsurv_artcd4dist], (4, 3, 2, 1))
  end
  if(size(par[:paedsurv_cd4dist])[1] != 53)
    par[:paedsurv_cd4dist] = permutedims(par[:paedsurv_cd4dist], (3, 2, 1))
  end
  global PROJ_YEARS = par[:ss][:PROJ_YEARS]
  global HIVSTEPS_PER_YEAR = par[:ss][:hiv_steps_per_year]

  global hAG_SPAN = par[:ss][:h_ag_span]
  global h_art_stage_dur = par[:ss][:h_art_stage_dur]
  global hAG_START = zeros(Int, hAG)
  for ha = 2:hAG
    hAG_START[ha] = hAG_START[ha-1] + hAG_SPAN[ha-1]
  end

  global DT = 1.0/HIVSTEPS_PER_YEAR
end

function simmodJ(par = nothing)::Dict
  if par === nothing
    par = prep_par()
  end
  prep_inputs!(par)
  out_dict = run_sim!(par)
  
  return out_dict
end
end
