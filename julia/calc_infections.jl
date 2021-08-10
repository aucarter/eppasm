function calc_infections_eppspectrum(pop, hivpop, artpop, r_ts, relinfectART, iota, incrr_sex, incrr_age,
  t_ART_start, DT, t, hts, hAG_START, hAG_SPAN, ts)

  # sum population sizes
  Xhivn_g = zeros(NG)
  Xhivn_incagerr = zeros(NG)
  Xhivp_noart = 0.0
  Xart = 0.0
  for g = 1:NG
    Xhivn_g[g] = 0.0
    Xhivn_incagerr[g] = 0.0
    for a = pIDX_15TO49:(pIDX_15TO49 + pAG_15TO49 - 1)
      Xhivn_g[g] += pop[t, HIVN, g, a]
      Xhivn_incagerr[g] += incrr_age[t, g, a] * pop[t, HIVN, g, a]
    end

    for ha = (hIDX_15TO49 + 1):hIDX_15TO49+hAG_15TO49+1

      # adjustment to first and last age group for partial year time step
      # calculation proportion of HIV population to include / exclude based on hivpop in single-year ages.
      if ha == hIDX_15TO49 
        hivp_ha = 0.0
        a = (hAG_START[ha] + 1):(hAG_START[ha]+hAG_SPAN[ha])
          hivp_ha += pop[t, HIVP, g, a]
        prop_include = (hivp_ha > 0) ? 1.0 - pop[t, HIVP, g, hAG_START[ha]] / hivp_ha * (1.0 - DT*hts) : 1.0
      elseif ha == hIDX_15TO49+hAG_15TO49
        hivp_ha = 0.0
        for a = (hAG_START[ha] + 1):(hAG_START[ha]+hAG_SPAN[ha])
          hivp_ha += pop[t, HIVP, g, a]
        end
        prop_include = (hivp_ha > 0) ? pop[t, HIVP, g, hAG_START[ha]] / hivp_ha * (1.0 - DT*hts) : 1.0
      else
        prop_include = 1.0
      end
      for hm = 1:hDS
        Xhivp_noart += hivpop[t, g, ha, hm] * prop_include
        if t >= t_ART_start 
          for hu = 1:hTS
            Xart += artpop[t, g, ha, hm, hu] * prop_include
          end
        end
      end
    end
  end # end loop over g
  Xhivn = Xhivn_g[MALE] + Xhivn_g[FEMALE]

  # adjust HIV negative population for partial year time step
  for g = 1:NG
    Xhivn -= pop[t, HIVN, g, pIDX_15TO49] * (1.0 - DT*hts)
    Xhivn += pop[t, HIVN, g, pIDX_15TO49 + pAG_15TO49 - 1] * (1.0 - DT*hts)
  end

  Xtot = Xhivn + Xhivp_noart + Xart
  prevcurr = (Xhivp_noart + Xart) / Xtot

  incrate15to49_ts[ts] = r_ts * (Xhivp_noart + relinfectART * Xart) / Xtot + iota


  # incidence by sex
  incrate15to49_g = zeros(NG)
  incrate15to49_g[MALE] = incrate15to49_ts[ts] * (Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + incrr_sex[t]*Xhivn_g[FEMALE])
  incrate15to49_g[FEMALE] = incrate15to49_ts[ts] * incrr_sex[t]*(Xhivn_g[MALE]+Xhivn_g[FEMALE]) / (Xhivn_g[MALE] + incrr_sex[t]*Xhivn_g[FEMALE])

  # annualized infections by age and sex
  infections_ts = zeros(NG, pAG)
  for g = 1:NG
    for a = 1:pAG
      infections_ts[g, a] = pop[t, HIVN, g, a] * incrate15to49_g[g] * incrr_age[t, g, a] * Xhivn_g[g] / Xhivn_incagerr[g]
    end
  end
  if isnan(prevcurr)
    error("Current prevalence is nan")
  end

  return infections_ts, prevcurr
end
