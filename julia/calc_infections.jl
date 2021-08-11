function calc_infections_eppspectrum!(mod_data, pop, hivpop, artpop, r_ts, iota, t, hts, ts)

  # sum population sizes
  Xhivn_g = zeros(NG)
  Xhivn_incagerr = zeros(NG)
  Xhivp_noart = 0.0
  Xart = 0.0
  for g = 1:NG
    Xhivn_g[g] = 0.0
    Xhivn_incagerr[g] = 0.0
    for a = pIDX_15TO49:(pIDX_15TO49 + pAG_15TO49 - 1)
      Xhivn_g[g] += pop[HIVN, g, a]
      Xhivn_incagerr[g] += mod_data[:incrr_age][a, g, t] * pop[HIVN, g, a]
    end

    for ha = (hIDX_15TO49 + 1):hIDX_15TO49+hAG_15TO49+1

      # adjustment to first and last age group for partial year time step
      # calculation proportion of HIV population to include / exclude based on hivpop in single-year ages.
      if ha == hIDX_15TO49 
        hivp_ha = 0.0
        a = (hAG_START[ha] + 1):(hAG_START[ha]+hAG_SPAN[ha])
          hivp_ha += pop[HIVP, g, a]
        prop_include = (hivp_ha > 0) ? 1.0 - pop[HIVP, g, hAG_START[ha]] / hivp_ha * (1.0 - DT*hts) : 1.0
      elseif ha == hIDX_15TO49+hAG_15TO49
        hivp_ha = 0.0
        for a = (hAG_START[ha] + 1):(hAG_START[ha]+hAG_SPAN[ha])
          hivp_ha += pop[HIVP, g, a]
        end
        prop_include = (hivp_ha > 0) ? pop[HIVP, g, hAG_START[ha]] / hivp_ha * (1.0 - DT*hts) : 1.0
      else
        prop_include = 1.0
      end
      for hm = 1:hDS
        Xhivp_noart += hivpop[g, ha, hm] * prop_include
        if t >= mod_data[:tARTstart]
          for hu = 1:hTS
            Xart += artpop[g, ha, hm, hu] * prop_include
          end
        end
      end
    end
  end # end loop over g
  Xhivn = Xhivn_g[MALE] + Xhivn_g[FEMALE]

  # adjust HIV negative population for partial year time step
  for g = 1:NG
    Xhivn -= pop[HIVN, g, pIDX_15TO49] * (1.0 - DT*hts)
    Xhivn += pop[HIVN, g, pIDX_15TO49 + pAG_15TO49 - 1] * (1.0 - DT*hts)
  end

  Xtot = Xhivn + Xhivp_noart + Xart
  prevcurr = (Xhivp_noart + Xart) / Xtot

  incrate15to49_ts[ts] = r_ts * (Xhivp_noart + mod_data[:relinfectART] * Xart) / Xtot + iota


  # incidence by sex
  incrate15to49_g = zeros(NG)
  incrate15to49_g[MALE] = incrate15to49_ts[ts] * (Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + mod_data[:incrr_sex][t]*Xhivn_g[FEMALE])
  incrate15to49_g[FEMALE] = incrate15to49_ts[ts] * mod_data[:incrr_sex][t]*(Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + mod_data[:incrr_sex][t]*Xhivn_g[FEMALE])

  # annualized infections by age and sex
  infections_ts = zeros(NG, pAG)
  for g = 1:NG
    for a = 1:pAG
      infections_ts[g, a] = pop[HIVN, g, a] * incrate15to49_g[g] * mod_data[:incrr_age][a, g, t] * Xhivn_g[g] / Xhivn_incagerr[g]
    end
  end
  if isnan(prevcurr)
    error("Current prevalence is nan")
  end

  return infections_ts, prevcurr
end
