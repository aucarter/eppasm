function hiv_mod_all(β, x, t, births_by_ha, out_dict)
  for hts = 0:(HIVSTEPS_PER_YEAR - 1)
    ts = (t - 2) * HIVSTEPS_PER_YEAR + hts + 1
    if β[:proj_steps][ts] >= β[:tsEpidemicStart]
      hiv_mod_one_step!(β, x, t, hts, ts, births_by_ha, out_dict)
    end
  end
end
  
function hiv_mod_one_step!(β, x, t, hts, ts, births_by_ha, out_dict)
  hivdeaths_ha = zeros(NG, hAG)
  disease_prog_mort!(β, x, hivdeaths_ha, t, out_dict)
  disease_transmission!(β, x, t, hts, ts, out_dict)
  if t >= β[:tARTstart]
    disease_treatment!(β, x, t, hts, births_by_ha, hivdeaths_ha, out_dict)
  end
  x[:hivpop] = x[:hivpop] + DT .* x[:grad]
  hivpop_ha = x[:pop][HIVP, :, :] * β[:age_binner]
  hivqx_ha = hivdeaths_ha ./ hivpop_ha
  replace!(hivqx_ha, NaN=>0.)
  hivqx_all_age = (hivqx_ha * β[:age_binner]')
  out_dict[:hivdeaths][t, :, :] = x[:pop][HIVP, :, :] .*  hivqx_all_age
  x[:pop][HIVP, :, :] .*= 1.0 .- hivqx_all_age
end

  
function disease_prog_mort!(β, x, hivdeaths_ha, t, out_dict)
  prog = β[:cd4_prog] .* x[:hivpop][:, :, 1:end - 1]
  x[:grad][:, :, 1:end - 1] -= prog
  x[:grad][:, :, 2:end] += prog 
  deaths = β[:cd4_mort] .* x[:hivpop]
  out_dict[:aidsdeaths_noart][t, :, :, :] .+= DT .* deaths
  hivdeaths_ha .+= DT .* dropdims(sum(deaths, dims = 3), dims = 3)
  x[:grad] .= x[:grad] - deaths
end

function disease_transmission!(β, x, t, hts, ts, out_dict)
  if eppmod != EPP_DIRECTINCID
    if eppmod == EPP_RSPLINE
      nothing
    else
      # Need to implement calc_rtrend_rt
      nothing
    end

    if incidmod == INCIDMOD_EPPSPEC 
      infections_ts, prevcurr = calc_infections_eppspectrum!(
        β, x, β[:rvec][ts], (β[:proj_steps][ts] == β[:tsEpidemicStart]) ? β[:iota] : 0.0,
        t, hts, ts, out_dict)
    end
    out_dict[:prev15to49_ts][ts] = prevcurr
    infections_a_dt = zeros(NG, pAG)
    infections_a_dt = DT .* infections_ts
    out_dict[:infections][t, :, :] .+= infections_a_dt
    x[:pop][HIVN, :, :] .-=  infections_a_dt
    x[:pop][HIVP, :, :] .+=  infections_a_dt
    infections_ha = infections_a_dt * β[:age_binner]
    x[:grad] .+= infections_ha .* β[:cd4_initdist]
  end
end

function calc_infections_eppspectrum!(β, x, r_ts, iota, t, hts, ts, out_dict)
  Xhivn_g = zeros(NG)
  Xhivn_incagerr = zeros(NG)
  Xhivp_noart = 0.0
  Xart = 0.0
  Xhivn_g .= dropdims(sum(x[:pop][HIVN, :, idx_15_49], dims = 2), dims = 2)
  Xhivn_incagerr .= dropdims(sum(β[:incrr_age][t, :, idx_15_49] .* x[:pop][HIVN, :, idx_15_49], dims = 2), dims = 2)
  Xhivp_noart = sum(x[:pop][HIVP, :, idx_15_49])
  Xart = sum(x[:artpop])
  # TODO: Implement this stuff for the edges
  # for g = 1:NG
  #   for ha = h_idx_15_49

  #     # adjustment to first and last age group for βtial year time step
  #     # calculation proportion of HIV x[:pop]ulation to include / exclude based on hivpop in single-year ages.
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
  #       if t >= β[:tARTstart]
  #         for hu = 1:hTS
  #           Xart += artpop[g, ha, hm, hu] * prop_include
  #         end
  #       end
  #     end
  #   end
  # end # end loop over g
  Xhivn = sum(Xhivn_g)

  # adjust HIV negative population for βtial year time step
  Xhivn -= sum(x[:pop][HIVN, :, pIDX_15TO49] .* (1.0 - DT*hts))
  Xhivn += sum(x[:pop][HIVN, :, pIDX_15TO49 + pAG_15TO49 - 1] .* (1.0 - DT*hts))
  Xtot = Xhivn + Xhivp_noart + Xart
  prevcurr = (Xhivp_noart + Xart) / Xtot

  out_dict[:incrate15to49_ts][ts] = r_ts * (Xhivp_noart + β[:relinfectART] * Xart) / Xtot + iota


  # incidence by sex
  incrate15to49_g = zeros(NG)
  incrate15to49_g[MALE] = out_dict[:incrate15to49_ts][ts] * (Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + β[:incrr_sex][t]*Xhivn_g[FEMALE])
  incrate15to49_g[FEMALE] = out_dict[:incrate15to49_ts][ts] * β[:incrr_sex][t]*(Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + β[:incrr_sex][t]*Xhivn_g[FEMALE])

  # annualized infections by age and sex
  infections_ts = zeros(NG, pAG)
  @. infections_ts = x[:pop][HIVN, :, :] * incrate15to49_g * β[:incrr_age][t, :, :] * Xhivn_g[:] / Xhivn_incagerr[:]
  if isnan(prevcurr)
    error("Current prevalence is nan")
  end

  # if sum(infections_ts .< -1000) > 0
  #   error("Really negative infections")
  # end

  return infections_ts, prevcurr
end
