function hiv_mod_all(p, x, t, births_by_ha, out_dict)
  for hts = 0:(p.HIVSTEPS_PER_YEAR - 1)
    ts = (t - 2) * p.HIVSTEPS_PER_YEAR + hts + 1
    if p.proj_steps[ts] >= p.tsEpidemicStart
      hiv_mod_one_step!(p, x, t, hts, ts, births_by_ha, out_dict)
    end
  end
end
  
function hiv_mod_one_step!(p, x, t, hts, ts, births_by_ha, out_dict)
  hivdeaths_ha = disease_prog_mort!(p, x, t, out_dict)
  disease_transmission!(p, x, t, hts, ts, out_dict)
  if t >= p.tARTstart
    disease_treatment!(p, x, t, hts, births_by_ha, hivdeaths_ha, out_dict)
  end
  x[:hivpop] = x[:hivpop] + p.DT .* x[:grad]
  hivpop_ha = x[:pop][p.HIVP, :, :] * p.age_binner
  hivqx_ha = hivdeaths_ha ./ hivpop_ha
  replace!(hivqx_ha, NaN=>0.)
  hivqx_all_age = (hivqx_ha * p.age_binner')
  out_dict[:hivdeaths][t, :, :] = x[:pop][p.HIVP, :, :] .*  hivqx_all_age
  x[:pop][p.HIVP, :, :] .*= 1.0 .- hivqx_all_age
end

  
function disease_prog_mort!(p, x, t, out_dict)
  x[:grad] .= 0.
  prog = p.cd4_prog .* x[:hivpop][:, :, 1:end - 1]
  x[:grad][:, :, 1:end - 1] -= prog
  x[:grad][:, :, 2:end] += prog 
  deaths = p.cd4_mort .* x[:hivpop]
  out_dict[:aidsdeaths_noart][t, :, :, :] .+= p.DT .* deaths
  hivdeaths_ha = p.DT .* dropdims(sum(deaths, dims = 3), dims = 3)
  x[:grad] .= x[:grad] - deaths
  return hivdeaths_ha
end

function disease_transmission!(p, x, t, hts, ts, out_dict)
  if p.eppmod != p.EPP_DIRECTINCID
    if p.eppmod == p.EPP_RSPLINE
      nothing
    else
      # Need to implement calc_rtrend_rt
      nothing
    end

    if p.incidmod == p.INCIDMOD_EPPSPEC 
      infections_ts, prevcurr = calc_infections_eppspectrum!(
        p, x, p.rvec[ts], (p.proj_steps[ts] == p.tsEpidemicStart) ? p.iota : 0.0,
        t, hts, ts, out_dict)
    end
    out_dict[:prev15to49_ts][ts] = prevcurr
    infections_a_dt = p.DT .* infections_ts
    out_dict[:infections][t, :, :] .+= infections_a_dt
    x[:pop][p.HIVN, :, :] .-=  infections_a_dt
    x[:pop][p.HIVP, :, :] .+=  infections_a_dt
    infections_ha = infections_ts * p.age_binner
    out_dict[:incid15to49][t] += p.DT * sum(infections_ha[p.h_idx_15_49]) 
    x[:grad] .+= infections_ha .* p.cd4_initdist
  end
end

function calc_infections_eppspectrum!(p, x, r_ts, iota, t, hts, ts, out_dict)
  Xhivn_g = zeros(p.NG)
  Xhivn_incagerr = zeros(p.NG)
  Xhivp_noart = 0.0
  Xart = 0.0
  Xhivn_g .= dropdims(sum(x[:pop][p.HIVN, :, p.idx_15_49], dims = 2), dims = 2)
  Xhivn_incagerr .= dropdims(sum(p.incrr_age[t, :, p.idx_15_49] .* x[:pop][p.HIVN, :, p.idx_15_49], dims = 2), dims = 2)
  prop_include = fill(1., p.NG, length(p.h_idx_15_49))
  prop_include[:, 1] .= 1.0 .- x[:pop][p.HIVP, :, p.hAG_START[(p.hIDX_15TO49 + 1)] + 1] ./  dropdims(sum(x[:pop][p.HIVP, :, (p.hAG_START[(p.hIDX_15TO49 + 1)] + 1):(p.hAG_START[(p.hIDX_15TO49 + 1)]+p.hAG_SPAN[(p.hIDX_15TO49 + 1)])], dims = 2), dims = 2) .* (1.0 - p.DT*hts)
  prop_include[:, end] .= x[:pop][p.HIVP, :, p.hAG_START[(p.hIDX_15TO49 + p.hAG_15TO49)]] ./ dropdims(sum(x[:pop][p.HIVP, :, (p.hAG_START[(p.hIDX_15TO49 + p.hAG_15TO49)]):(p.hAG_START[(p.hIDX_15TO49 + p.hAG_15TO49) + 1]+p.hAG_SPAN[(p.hIDX_15TO49 + p.hAG_15TO49) + 1])], dims = 2), dims = 2) .* (1.0 - p.DT*hts)
  replace!(prop_include, NaN=>0.)
  Xhivp_noart = sum(prop_include .* x[:hivpop][:,p.h_idx_15_49, :])
  Xart = sum(prop_include .* x[:artpop])
  Xhivn = sum(Xhivn_g)

  # adjust HIV negative population for ptial year time step
  Xhivn -= sum(x[:pop][p.HIVN, :, p.pIDX_15TO49] .* (1.0 - p.DT*hts))
  Xhivn += sum(x[:pop][p.HIVN, :, p.pIDX_15TO49 + p.pAG_15TO49 - 1] .* (1.0 - p.DT*hts))
  Xtot = Xhivn + Xhivp_noart + Xart
  prevcurr = (Xhivp_noart + Xart) / Xtot

  out_dict[:incrate15to49_ts][ts] = r_ts * (Xhivp_noart + p.relinfectART * Xart) / Xtot + iota


  # incidence by sex
  incrate15to49_g = zeros(p.NG)
  incrate15to49_g[p.MALE] = out_dict[:incrate15to49_ts][ts] * (Xhivn_g[p.MALE]+Xhivn_g[p.FEMALE]) / (Xhivn_g[p.MALE] + p.incrr_sex[t]*Xhivn_g[p.FEMALE])
  incrate15to49_g[p.FEMALE] = out_dict[:incrate15to49_ts][ts] * p.incrr_sex[t]*(Xhivn_g[p.MALE]+Xhivn_g[p.FEMALE]) / (Xhivn_g[p.MALE] + p.incrr_sex[t]*Xhivn_g[p.FEMALE])

  # annualized infections by age and sex
  infections_ts = zeros(p.NG, p.pAG)
  @. infections_ts = x[:pop][p.HIVN, :, :] * incrate15to49_g * p.incrr_age[t, :, :] * Xhivn_g[:] / Xhivn_incagerr[:]
  @assert !isnan(prevcurr)

  return infections_ts, prevcurr
end
