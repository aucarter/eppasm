function hiv_mod_all(par, t, anyelig_idx, cd4elig_idx, artpop, hivpop, aidsdeaths_noart, everARTelig_idx, births_by_ha, pop, age_binner)
  for hts = 0:(HIVSTEPS_PER_YEAR - 1)
    ts = (t - 2) * HIVSTEPS_PER_YEAR + hts + 1
    if par[:proj_steps][ts] >= par[:tsEpidemicStart]
      hiv_mod_one_step!(par, t, hts, ts,  anyelig_idx, cd4elig_idx, artpop, hivpop, aidsdeaths_noart, everARTelig_idx, births_by_ha, pop, age_binner)
    end
  end
end
  
function hiv_mod_one_step!(par, t, hts, ts, anyelig_idx, cd4elig_idx, artpop, hivpop, aidsdeaths_noart, everARTelig_idx, births_by_ha, pop, age_binner)
  hivdeaths_ha = zeros(NG, hAG)
  grad = zeros(NG, hAG, hDS)
  
  disease_prog_mort!(par, hivpop, hivdeaths_ha, grad, aidsdeaths_noart, t, everARTelig_idx)
  disease_transmission!(par, pop, hivpop, artpop, grad, t, hts, ts, age_binner)
  if t >= par[:tARTstart]
    disease_treatment!(
      par, pop, hivpop, artpop, grad, t, hts, births_by_ha, hivdeaths_ha, everARTelig_idx,
      anyelig_idx, cd4elig_idx, age_binner
    )
  end
  hivpop = hivpop + DT .* grad
  hivpop_ha = pop[HIVP, :, :] * age_binner
  hivqx_ha = hivdeaths_ha ./ hivpop_ha
  replace!(hivqx_ha, NaN=>0.)
  hivqx_all_age = (hivqx_ha * age_binner')
  hivdeaths[t, :, :] = pop[HIVP, :, :] .*  hivqx_all_age
  pop[HIVP, :, :] .*= 1.0 .- hivqx_all_age
end

  
function disease_prog_mort!(par, hivpop, hivdeaths_ha, grad, aidsdeaths_noart, t, everARTelig_idx)
  prog = par[:cd4_prog] .* hivpop[:, :, 1:end - 1]
  grad[:, :, 1:end - 1] -= prog
  grad[:, :, 2:end] += prog 
  deaths = par[:cd4_mort] .* hivpop
  aidsdeaths_noart[t, :, :, :] .+= DT .* deaths
  hivdeaths_ha .+= DT .* dropdims(sum(deaths, dims = 3), dims = 3)
  grad = grad - deaths
end

function disease_transmission!(par, pop, hivpop, artpop, grad, t, hts, ts, age_binner)
  if eppmod != EPP_DIRECTINCID
    if eppmod == EPP_RSPLINE
      nothing
    else
      # Need to implement calc_rtrend_rt
      nothing
    end

    if incidmod == INCIDMOD_EPPSPEC 
      infections_ts, prevcurr = calc_infections_eppspectrum!(
        par, pop, hivpop, artpop, rvec[ts], (par[:proj_steps][ts] == par[:tsEpidemicStart]) ? par[:iota] : 0.0,
        t, hts, ts)
    end
    global prev15to49_ts[ts] = prevcurr
    infections_a_dt = zeros(NG, pAG)
    infections_a_dt = DT .* infections_ts
    pop[HIVN, :, :] .-=  infections_a_dt
    pop[HIVP, :, :] .+=  infections_a_dt
    infections_ha = infections_a_dt * age_binner
    grad .+= infections_ha .* par[:cd4_initdist]
  end
end

function calc_infections_eppspectrum!(par, pop, hivpop, artpop, r_ts, iota, t, hts, ts)
  Xhivn_g = zeros(NG)
  Xhivn_incagerr = zeros(NG)
  Xhivp_noart = 0.0
  Xart = 0.0
  Xhivn_g .= dropdims(sum(pop[HIVN, :, idx_15_49], dims = 2), dims = 2)
  Xhivn_incagerr .= dropdims(
    sum(
      par[:incrr_age][t, :, idx_15_49] .* pop[HIVN, :, idx_15_49], dims = 2
    ), dims = 2
  )
  Xhivp_noart = sum(pop[HIVP, :, idx_15_49])
  Xart = sum(artpop)
  # TODO: Implement this stuff for the edges
  # for g = 1:NG
  #   for ha = h_idx_15_49

  #     # adjustment to first and last age group for partial year time step
  #     # calculation proportion of HIV population to include / exclude based on hivpop in single-year ages.
  #     if ha == (hIDX_15TO49 + 1)
  #       hivp_ha = 0.0
  #       for a = (hAG_START[ha] + 1):(hAG_START[ha]+hAG_SPAN[ha])
  #         hivp_ha += pop[HIVP, g, a]
  #       end
  #       prop_include = (hivp_ha > 0) ? 1.0 - pop[HIVP, g, hAG_START[ha] + 1] / hivp_ha * (1.0 - DT*hts) : 1.0
  #     elseif ha == (hIDX_15TO49 + hAG_15TO49)
  #       hivp_ha = 0.0
  #       for a = (hAG_START[ha] + 1):(hAG_START[ha]+hAG_SPAN[ha])
  #         hivp_ha += pop[HIVP, g, a]
  #       end
  #       prop_include = (hivp_ha > 0) ? pop[HIVP, g, hAG_START[ha] + 1] / hivp_ha * (1.0 - DT*hts) : 1.0
  #     else
  #       prop_include = 1.0
  #     end
  #     for hm = 1:hDS
  #       Xhivp_noart += hivpop[g, ha, hm] * prop_include
  #       if t >= par[:tARTstart]
  #         for hu = 1:hTS
  #           Xart += artpop[g, ha, hm, hu] * prop_include
  #         end
  #       end
  #     end
  #   end
  # end # end loop over g
  Xhivn = sum(Xhivn_g)

  # adjust HIV negative population for partial year time step
  Xhivn -= sum(pop[HIVN, :, pIDX_15TO49] .* (1.0 - DT*hts))
  Xhivn += sum(pop[HIVN, :, pIDX_15TO49 + pAG_15TO49 - 1] .* (1.0 - DT*hts))
  Xtot = Xhivn + Xhivp_noart + Xart
  prevcurr = (Xhivp_noart + Xart) / Xtot

  incrate15to49_ts[ts] = r_ts * (Xhivp_noart + par[:relinfectART] * Xart) / Xtot + iota


  # incidence by sex
  incrate15to49_g = zeros(NG)
  incrate15to49_g[MALE] = incrate15to49_ts[ts] * (Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + par[:incrr_sex][t]*Xhivn_g[FEMALE])
  incrate15to49_g[FEMALE] = incrate15to49_ts[ts] * par[:incrr_sex][t]*(Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + par[:incrr_sex][t]*Xhivn_g[FEMALE])

  # annualized infections by age and sex
  infections_ts = zeros(NG, pAG)
  @. infections_ts = pop[HIVN, :, :] * incrate15to49_g * par[:incrr_age][t, :, :] * Xhivn_g[:] / Xhivn_incagerr[:]
  if isnan(prevcurr)
    error("Current prevalence is nan")
  end

  return infections_ts, prevcurr
end
