using LinearAlgebra
using RCall
using Serialization

include("pop_project.jl")
include("hiv_mod.jl")
include("calc_infections.jl")

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
const hIDX_CD4_350 = 2
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

function prep_mod_data()
  # mod_data = rcopy(R"""
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
  mod_data  = deserialize("test")
  return mod_data
end

function sim_one_year!(
  t, mod_data, pop, pop_ts, hivpop, hivpop_ts, artpop, artpop_ts,
  everARTelig_idx, age_binner, max_age_binner)

  # Project forward the population
  births, births_by_ha, last_hivpop, last_artpop = pop_project_one_step!(mod_data, t, pop, hivpop, artpop,  age_binner, max_age_binner)
  
  ## HIV model simulation
  cd4elig_idx = artcd4elig_idx[t] 
  anyelig_idx = (specpop_percelig[t] > 0 || pw_artelig[t] > 0) ? 1 : (who34percelig > 0) ? hIDX_CD4_350 : cd4elig_idx
  everARTelig_idx = (anyelig_idx < everARTelig_idx) ? anyelig_idx : everARTelig_idx
  if mod_data[:proj_steps][(t - 2) * HIVSTEPS_PER_YEAR + 1] >= mod_data[:tsEpidemicStart]
    for hts = 0:(HIVSTEPS_PER_YEAR - 1)
      ts = (t - 2) * HIVSTEPS_PER_YEAR + hts + 1
      if mod_data[:proj_steps][ts] >= mod_data[:tsEpidemicStart]
        hiv_mod_one_step!(mod_data, t, hts, anyelig_idx, cd4elig_idx, artpop, hivpop, aidsdeaths_noart, everARTelig_idx, births_by_ha, pop, age_binner)
      end
    end
  end
  
  # Adjust to a target population
  if mod_data[:popadjust]
    adjust_pop!(mod_data, t, pop, hivpop, artpop)
  end

  # Prevalence in pregnant women
  preg_women_prev_one_step!(mod_data, t, pregprev, pop_ts, pop, hivpop, artpop, births, last_hivpop, last_artpop)

  # Update output timeseries objects
  pop_ts[t, :, :, :] = pop
  hivpop_ts[t, :, :, :] = hivpop
  artpop_ts[t, :, :, :, :] = artpop
end

function run_sim!(mod_data)
  # Setup baseline
  pop = zeros(pDS, NG, pAG)
  pop[HIVN, :, :] = mod_data[:basepop]'
  pop[HIVP, :, :] .= 0.
  hivpop = zeros(NG, hAG, hDS)
  artpop = zeros(NG, hAG, hDS, hTS)
  if(size(mod_data[:targetpop])[1] == 66)
    mod_data[:targetpop] = permutedims(mod_data[:targetpop], [3, 2, 1])
  end
  # Track updates to these state objects in each year
  pop_ts = zeros(PROJ_YEARS, pDS, NG, pAG)
  hivpop_ts = zeros(PROJ_YEARS, NG, hAG, hDS)
  artpop_ts = zeros(PROJ_YEARS, NG, hAG, hDS, hTS)
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

  # Simulate all years
  for t = 2:mod_data[:SIM_YEARS]
    sim_one_year!(t, mod_data, pop, pop_ts, hivpop, hivpop_ts, artpop, artpop_ts,
      everARTelig_idx, age_binner, max_age_binner)
  end

  # Calculate outputs
  hivn15to49 = dropdims(sum(pop_ts[:, HIVN, :, pIDX_15TO49:(pIDX_15TO49 + pAG_15TO49 - 1)], dims = (2, 3, 4)), dims = (2, 3))
  hivp15to49 = dropdims(sum(pop_ts[:, HIVP, :, pIDX_15TO49:(pIDX_15TO49 + pAG_15TO49 - 1)], dims = (2, 3, 4)), dims = (2, 3))
  prev15to49 = hivp15to49 ./ (hivn15to49 + hivp15to49)
  incid15to49[2:end] ./= hivn15to49[1:end - 1]

  return pop_ts, hivpop_ts, artpop_ts, prev15to49
end

function simmodJ(mod_data = nothing)
  if mod_data === nothing
    mod_data = prep_mod_data()
  end
  if(size(mod_data[:Sx])[1] == 66)
    mod_data[:Sx] = permutedims(mod_data[:Sx], [3, 2, 1])
  end
  if(size(mod_data[:netmigr])[1] == 66)
    mod_data[:netmigr] = permutedims(mod_data[:netmigr], [3, 2, 1])
  end
  if(size(mod_data[:cd4_mort])[1] != 2)
    mod_data[:cd4_mort] = permutedims(mod_data[:cd4_mort], (3, 2, 1))
  end
  if(size(mod_data[:cd4_prog])[1] != 2)
    mod_data[:cd4_prog] = permutedims(mod_data[:cd4_prog], (3, 2, 1))
  end
  global PROJ_YEARS = mod_data[:ss][:PROJ_YEARS]
  global HIVSTEPS_PER_YEAR = mod_data[:ss][:hiv_steps_per_year]

  global hAG_SPAN = mod_data[:ss][:h_ag_span]
  global h_art_stage_dur = mod_data[:ss][:h_art_stage_dur]
  global hAG_START = zeros(Int, hAG)
  for ha = 2:hAG
    hAG_START[ha] = hAG_START[ha-1] + hAG_SPAN[ha-1]
  end
  global entrantprev = mod_data[:entrantprev]'
  global use_entrantprev = true
  global verttrans_lag = mod_data[:verttrans_lag]
  global paedsurv_lag = mod_data[:paedsurv_lag]
  global artcd4elig_idx = mod_data[:artcd4elig_idx]
  global specpop_percelig = mod_data[:specpop_percelig]
  global pw_artelig = mod_data[:pw_artelig]
  global who34percelig = mod_data[:who34percelig]

  global prev15to49_ts = zeros((PROJ_YEARS-1) * HIVSTEPS_PER_YEAR)
  global incrate15to49_ts = zeros((PROJ_YEARS-1) * HIVSTEPS_PER_YEAR)
  global infections = zeros(PROJ_YEARS, NG, pAG)
  global hivdeaths = zeros(PROJ_YEARS, NG, pAG)
  global natdeaths = zeros(PROJ_YEARS, NG, pAG)
  global incid15to49 = zeros(PROJ_YEARS)
  global aidsdeaths_noart = zeros(PROJ_YEARS, NG, hAG, hDS)
  global aidsdeaths_art = zeros(PROJ_YEARS, NG, hAG, hDS, hTS)
  global artinit = zeros(PROJ_YEARS, NG, hAG, hDS)
  global pregprevlag = zeros(PROJ_YEARS)
  global pregprev = zeros(PROJ_YEARS)
  global popadjust = zeros(PROJ_YEARS, NG, pAG)
  global entrantprev_out = zeros(PROJ_YEARS)

  global DT = 1.0/HIVSTEPS_PER_YEAR
  if eppmod == EPP_RSPLINE
    global rvec = mod_data[:rvec]
  end

  pop_ts, hivpop_ts, artpop_ts, prev15to49 = run_sim!(mod_data)
  
  out_dict = Dict(
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
    :rvec => rvec,
    :prev15to49 => prev15to49,
    :pregprev => pregprev,
    :incid15to49 => incid15to49,
    :entrantprev_out => entrantprev_out,
    :pop => pop_ts
  )
  return out_dict
end

mod_data = prep_mod_data()
@time mod_dict = simmodJ(mod_data)
mod_data = prep_mod_data()
@time mod_dict = simmodJ(mod_data)
using Plots
plot_pop_ta = sum(mod_dict[:pop], dims = (2, 3))[:, 1, 1, :]
heatmap(plot_pop_ta')