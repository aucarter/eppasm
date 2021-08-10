function pop_project_one_step(t, pop, hivpop, artpop)
  ## Population projection
  # Age pop one year
  pop[t, :, :, 2:end] = pop[t - 1, :, :, 1:(end - 1)]
  pop[t, :, :, end] = pop[t - 1, :, :, end]

  # Calculate probability of aging from one HIV positive age bin to the next
  if sum(pop[t-1, HIVP, :, :]) > 0
    hiv_ag_prob = (pop[t-1, HIVP, :, :] * max_age_binner) ./ (pop[t-1, HIVP, :, :] * age_binner)
    replace!(hiv_ag_prob, NaN=>0.)
    hiv_ag_prob[:, end] .= 0
  else
    hiv_ag_prob = zeros(NG, hAG)
  end

  # Age the HIV positive population
  @. hivpop[t, :, 2:end, :] =  (1-hiv_ag_prob[:, 2:end]) * hivpop[t-1, :, 2:end, : ] + hiv_ag_prob[:, 1:end-1] * hivpop[t-1, :, 1:end-1, :]
  if t > t_ART_start
    @. artpop[t, :, 2:end, :, :] =  (1-hiv_ag_prob[:, 2:end]) * artpop[t-1, :, 2:end, :, :] + hiv_ag_prob[:, 1:end-1] * artpop[t-1, :, 1:end-1, :, :]
  end

  # Lagged birth to youngest age group
  for g = 1:NG
    if use_entrantprev
      entrant_prev = entrantprev[t, g]
    else
      entrant_prev = pregprevlag[t-1] * verttrans_lag[t-1] * paedsurv_lag[t-1]
    end

    if bin_popadjust
      pop[t, HIVN, g, 1] =  entrantpop[t-1, g] * (1.0-entrant_prev)
      paedsurv_g = entrantpop[t-1, g] * entrant_prev
    else
      pop[t, HIVN, g, 1] = birthslag[t-1, g] * cumsurv[t-1, g] * (1.0-entrant_prev / paedsurv_lag[t-1]) + cumnetmigr[t-1, g] * (1.0-pregprevlag[t-1] * netmig_hivprob)
      paedsurv_g = birthslag[t-1, g] * cumsurv[t-1, g] * entrant_prev + cumnetmigr[t-1, g] * entrant_prev    
    end

    pop[t, HIVP, g, 1] = paedsurv_g
    entrantprev_out[t] = (pop[t, HIVP, MALE, 1] + pop[t, HIVP, FEMALE, 1]) / (pop[t, HIVN, MALE, 1] + pop[t, HIVN, FEMALE, 1] + pop[t, HIVP, MALE, 1] + pop[t, HIVP, FEMALE, 1])

    for hm = 1:hDS
      hivpop[t, g, 1, hm] = (1-hiv_ag_prob[g, 1]) * hivpop[t-1, g, 1, hm] + paedsurv_g * paedsurv_cd4dist[t, g, hm] * (1.0 - entrantartcov[t, g])
      if t > t_ART_start
        for hu = 1:hTS
          artpop[t, g, 1, hm, hu] = (1-hiv_ag_prob[g, 1]) * artpop[t-1, g, 1, hm, hu]
          artpop[t, g, 1, hm, hu] += paedsurv_g * paedsurv_artcd4dist[t, g, hm, hu] * entrantartcov[t, g]
        end
      end
    end
  end

  # non-HIV mortality and migration
  hivpop_ha = pop[t, HIVP, :, :] * age_binner
  qx = 1 .- Sx[t, :, :]
  ndeaths = pop[t, HIVN, :, :] .* qx
  pop[t, HIVN, :, :] .-= ndeaths
  hdeaths = pop[t, HIVP, :, :] .* qx
  pop[t, HIVP, :, :] .-= hdeaths

  # net migration
  migrate = @. netmigr[t, :, :] * (1+Sx[t, :, :])/2.0 / (pop[t, HIVN, :, :] + pop[t, HIVP, :, :])
  pop[t, HIVN, :, :] .*= 1 .+ migrate
  hmig = migrate .* pop[t, HIVP, :, :]
  deathsmig = hmig * age_binner
  pop[t, HIVP, :, :] .+= hmig
  deathmigrate = deathsmig ./ hivpop_ha
  replace!(deathmigrate, NaN=>0.)
  hivpop[t, :, :, :] .*= 1 .+ deathmigrate
  if(t > t_ART_start)
    artpop[t, :, :, :, :] .*= 1 .+ deathmigrate
  end
  natdeaths[t, :, :] .= ndeaths .+ hdeaths

  # Fertility
  global births = 0.0
  global births_by_ha = zeros(hAG_FERT)
  for m = 1:pDS
    a = pIDX_FERT
    for ha = (hIDX_FERT + 1):(hIDX_FERT+hAG_FERT)
      for i = 1:hAG_SPAN[ha]
        births_by_ha[ha-hIDX_FERT] += (pop[t-1, m, FEMALE, a] + pop[t, m, FEMALE, a])/2 * asfr[t, a]
        a += 1
      end
    end
  end
  for ha = (hIDX_FERT + 1):hAG_FERT
    births += births_by_ha[ha-hIDX_FERT]
  end
  if t + AGE_START < PROJ_YEARS 
    for g = 1:NG
      birthslag[(t + AGE_START-1), g] = srb[t, g] * births
    end
  end

  return births
end

function preg_women_prev_one_step(t, pregprev, pop, artpop, births)
  # Prevalence among pregnant women
  hivbirths = 0.
  for ha = (hIDX_FERT + 1):(hIDX_FERT + hAG_FERT)
    hivn_ha = 0.
    frr_hivpop_ha = 0.
    for a = (hAG_START[ha] + 1):(hAG_START[ha]+hAG_SPAN[ha])
      hivn_ha += (pop[t-1, HIVN, FEMALE, a] + pop[t, HIVN, FEMALE, a])/2
    end
    for hm = 1:hDS
      frr_hivpop_ha += frr_cd4[t, ha-hIDX_FERT, hm] * (hivpop[t-1, FEMALE, ha, hm]+hivpop[t, FEMALE, ha, hm])/2
      if t == t_ART_start 
        for hu = 1:hTS
          frr_hivpop_ha += frr_art[t, ha-hIDX_FERT, hm, hu] * artpop[t, FEMALE, ha, hm, hu]/2
        end
      elseif t > t_ART_start
        for hu = 1:hTS
          frr_hivpop_ha += frr_art[t, ha-hIDX_FERT, hm, hu] * (artpop[t-1, FEMALE, ha, hm, hu]+artpop[t, FEMALE, ha, hm, hu])/2
        end
      end
    end
  end

  pregprev[t] = hivbirths/births
  if (t + AGE_START) < PROJ_YEARS
    pregprevlag[t + AGE_START] = pregprev[t]
  end
end