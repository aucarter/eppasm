using LinearAlgebra
using RCall


include("pop_project.jl")
include("hiv_mod.jl")
include("calc_infections.jl")

const AGE_START = 15
const NG = 2 # Genders
const pAG = 66 # Ages
const pDS = 2 # Disease states
const pIDX_FERT = 1
const pAG_FERT = 35
const pIDX_15TO49 = 1
const pAG_15TO49 = 35
const pIDX_15PLUS = 1
const pAG_15PLUS = 66
const hAG = 9
const hDS = 7
const hTS = 3
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

function run_sim!(pop, hivpop, artpop, everARTelig_idx)
  for t = 2:SIM_YEARS
    pop_project_one_step(t, pop, hivpop, artpop)
    ## HIV model simulation
    cd4elig_idx = artcd4elig_idx[t] 
    anyelig_idx = (specpop_percelig[t] > 0 || pw_artelig[t] > 0) ? 1 : (who34percelig > 0) ? hIDX_CD4_350 : cd4elig_idx
    everARTelig_idx = (anyelig_idx < everARTelig_idx) ? anyelig_idx : everARTelig_idx
    
    for hts = 0:(HIVSTEPS_PER_YEAR - 1)
      hiv_mod_one_step(t, hts, anyelig_idx, cd4elig_idx, artpop, hivpop, aidsdeaths_noart, everARTelig_idx)
    end
    
    if bin_popadjust
      for g = 1:NG
        a = 1
        for ha = 1:hAG
          popadj_ha = 0.
          hivpop_ha = 0.
          for i = 1:hAG_SPAN[ha]
            hivpop_ha += pop[t, HIVP, g, a]
            popadjrate_a = popadjust[t, g, a] = targetpop[t, g, a] / (pop[t, HIVN, g, a] + pop[t, HIVP, g, a])
            pop[t, HIVN, g, a] *= popadjrate_a
            hpopadj_a = (popadjrate_a-1.0) * pop[t, HIVP, g, a]
            popadj_ha += hpopadj_a
            pop[t, HIVP, g, a] += hpopadj_a
            a += 1
          end
          popadjrate_ha = hivpop_ha > 0 ? popadj_ha / hivpop_ha : 0.0
          for hm = 1:hDS
            hivpop[t, g, ha, hm] *= 1+popadjrate_ha
            if t >= t_ART_start
              for hu = 1:hTS
                artpop[t, g, ha, hm, hu] *= 1+popadjrate_ha
              end
            end
          end
        end
      end
    end

    preg_women_prev_one_step(t, pregprev, pop, artpop, births)

    hivn15to49[t] = sum(pop[t, HIVN, :, pIDX_15TO49:(pIDX_15TO49 + pAG_15TO49 - 1)])
    hivp15to49[t] = sum(pop[t, HIVP, :, pIDX_15TO49:(pIDX_15TO49 + pAG_15TO49 - 1)])
    prev15to49[t] = hivp15to49[t]/(hivn15to49[t] + hivp15to49[t])
    incid15to49[t] /= hivn15to49[t-1]
  end
  hivn15to49 = dropdims(sum(pop[:, HIVN, :, pIDX_15TO49:(pIDX_15TO49 + pAG_15TO49 - 1)], dims = (2, 3, 4)), dims = (2, 3, 4))
  # hivp15to49[t] = sum(pop[t, HIVP, :, pIDX_15TO49:(pIDX_15TO49 + pAG_15TO49 - 1)])
  # prev15to49[t] = hivp15to49[t]/(hivn15to49[t] + hivp15to49[t])
  # incid15to49[t] /= hivn15to49[t-1]
end

function simmodJ(mod_data = nothing)
  if mod_data === nothing
    mod_data = rcopy(R"""
    library(eppasm)
    pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
    bw <- prepare_spec_fit(pjnz, proj.end=2022.5)
    bw_fp <- attr(bw$Urban, "specfp")
    bw_fp <- prepare_rhybrid(bw_fp)
    bw_fp$ancsitedata <- TRUE
    bw_fp$ancrt <- "census"
    bw_fp$ancrtsite.beta <- 0
    bw_fp$logitiota <- TRUE
    bw_theta <- c(-0.407503322169364, -2.76794181367538, -1.26018073624346, 1995.96447776502,
                  -0.00307437171215574, 0.0114118307148102, 0.00760958379603691, 0.02,
                  2.24103194827232, -0.0792123921862689, -5.01917961803606, 0.359444135205712,
                  -6.10051517060137)
    param <- fnCreateParam(bw_theta, bw_fp)
    bw_fp <- update(bw_fp, list=param)
    bw_fp$eppmodInt <- match(bw_fp$eppmod, c("rtrend", "directincid"), nomatch=0) # 0: r-spline;
    bw_fp$incidmodInt <- match(bw_fp$incidmod, c("eppspectrum"))-1L
    bw_fp
    """)
  end
  global basepop = mod_data[:basepop]'
  global Sx = permutedims(mod_data[:Sx], [3, 2, 1])
  global netmigr = permutedims(mod_data[:netmigr], [3, 2, 1])
  global asfr = mod_data[:asfr]'
  global srb = mod_data[:srb]'
  global birthslag = mod_data[:birthslag]'
  global f_idx = mod_data[:ss][:p_fert_idx]
  global PROJ_YEARS = mod_data[:ss][:PROJ_YEARS]
  global HIVSTEPS_PER_YEAR = mod_data[:ss][:hiv_steps_per_year]

  global hAG_SPAN = mod_data[:ss][:h_ag_span]
  global h_art_stage_dur = mod_data[:ss][:h_art_stage_dur]
  global hAG_START = zeros(Int, hAG)
  for ha = 2:hAG
    hAG_START[ha] = hAG_START[ha-1] + hAG_SPAN[ha-1]
  end
  global SIM_YEARS = mod_data[:SIM_YEARS]
  global entrantprev = mod_data[:entrantprev]'
  global use_entrantprev = true
  global verttrans_lag = mod_data[:verttrans_lag]
  global paedsurv_lag = mod_data[:paedsurv_lag]
  global bin_popadjust = mod_data[:popadjust]
  global targetpop = permutedims(mod_data[:targetpop], [3, 2, 1])
  global entrantpop = mod_data[:entrantpop]'
  global paedsurv_cd4dist = permutedims(mod_data[:paedsurv_cd4dist], [3, 2, 1])
  global paedsurv_artcd4dist = permutedims(mod_data[:paedsurv_artcd4dist], [4, 3, 2, 1])
  global entrantartcov = mod_data[:entrantartcov]'
  global rspline_rvec = mod_data[:rvec]
  global projsteps = mod_data[:proj_steps]
  global tsEpidemicStart = mod_data[:tsEpidemicStart]
  global relinfectART = mod_data[:relinfectART]
  global iota = mod_data[:iota]
  global incrr_sex = mod_data[:incrr_sex]
  global incrr_age = permutedims(mod_data[:incrr_age], (3, 2, 1))
  global t_ART_start = mod_data[:tARTstart]
  global artnum15plus = mod_data[:art15plus_num]'
  global art15plus_isperc = mod_data[:art15plus_isperc]'
  global artcd4elig_idx = mod_data[:artcd4elig_idx]
  global specpop_percelig = mod_data[:specpop_percelig]
  global pw_artelig = mod_data[:pw_artelig]
  global who34percelig = mod_data[:who34percelig]
  global art_dropout = mod_data[:art_dropout]
  global median_cd4init = mod_data[:median_cd4init]
  global med_cd4init_cat = mod_data[:med_cd4init_cat]
  global med_cd4init_input = mod_data[:med_cd4init_input]
  global art_alloc_method = mod_data[:art_alloc_method]
  global art_alloc_mxweight = mod_data[:art_alloc_mxweight][1]
  global scale_cd4_mort = mod_data[:scale_cd4_mort]
  global cd4_initdist = permutedims(mod_data[:cd4_initdist], [3, 2, 1])
  global cd4_prog = permutedims(mod_data[:cd4_prog], [3, 2, 1])
  global cd4_mort = permutedims(mod_data[:cd4_mort], [3, 2, 1])
  global art_mort = permutedims(mod_data[:art_mort], [4, 3, 2, 1])
  global artmx_timerr = mod_data[:artmx_timerr]'
  global frr_cd4 = permutedims(mod_data[:frr_cd4], [3, 2, 1])
  global frr_art = permutedims(mod_data[:frr_art], [4, 3, 2, 1])

  global prev15to49_ts = zeros((PROJ_YEARS-1) * HIVSTEPS_PER_YEAR)
  global incrate15to49_ts = zeros((PROJ_YEARS-1) * HIVSTEPS_PER_YEAR)
  global prev15to49 = zeros(PROJ_YEARS)
  global infections = zeros(PROJ_YEARS, NG, pAG)
  global hivdeaths = zeros(PROJ_YEARS, NG, pAG)
  global natdeaths = zeros(PROJ_YEARS, NG, pAG)
  global incid15to49 = zeros(PROJ_YEARS)
  global aidsdeaths_noart = zeros(PROJ_YEARS, NG, hAG, hDS)
  global aidsdeaths_art = zeros(PROJ_YEARS, NG, hAG, hDS, hTS)
  global artinit = zeros(PROJ_YEARS, NG, hAG, hDS)
  global pop = zeros(PROJ_YEARS, pDS, NG, pAG)
  global hivn15to49 = zeros(PROJ_YEARS)
  global hivp15to49 = zeros(PROJ_YEARS)
  global hivpop = zeros(PROJ_YEARS, NG, hAG, hDS)
  global artpop = zeros(PROJ_YEARS, NG, hAG, hDS, hTS)
  global pregprevlag = zeros(PROJ_YEARS)
  global pregprev = zeros(PROJ_YEARS)
  global popadjust = zeros(PROJ_YEARS, NG, pAG)
  global entrantprev_out = zeros(PROJ_YEARS)

  everARTelig_idx = hDS
  global DT = 1.0/HIVSTEPS_PER_YEAR
  
  global age_binner = zeros(sum(hAG_SPAN), length(hAG_SPAN))
  global max_age_binner = zeros(sum(hAG_SPAN), length(hAG_SPAN))
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
  
  pop[1, HIVN, :, :] = basepop
  pop[1, HIVP, :, :] .= 0.
  hivn15to49[1] = sum(basepop[:, pIDX_15TO49:(pIDX_15TO49 + pAG_15TO49 - 1)])

  run_sim!(pop, hivpop, artpop, everARTelig_idx)

  # return hivpop, artpop, infections, hivdeaths, natdeaths, aidsdeaths_noart, aidsdeaths_art, popadjust, artinit, pregprevlag, incrate15to49_ts, prev15to49_ts, rvec, prev15to49, pregprev, incid15to49, entrantprev_out

end
@time simmodJ()
# hivpop, artpop, infections, hivdeaths, natdeaths, aidsdeaths_noart, aidsdeaths_art, popadjust, artinit, pregprevlag, incrate15to49_ts, prev15to49_ts, rvec, prev15to49, pregprev, incid15to49, entrantprev_out = simmodJ()
# plot_pop_ta = sum(pop, dims = (2, 3))[:, 1, 1, :]
# heatmap(plot_pop_ta')