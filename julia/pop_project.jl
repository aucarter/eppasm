function pop_project_one_step!(par, t, pop, hivpop, artpop, age_binner, max_age_binner)
  ## Population projection
  # Calculate probability of aging from one HIV positive age bin to the next
  if sum(pop[HIVP, :, :]) > 0
    hiv_ag_prob = (pop[HIVP, :, :] * max_age_binner) ./ (pop[HIVP, :, :] * age_binner)
    replace!(hiv_ag_prob, NaN=>0.)
    hiv_ag_prob[:, end] .= 0
  else
    hiv_ag_prob = zeros(NG, hAG)
  end

  # Save last year fertile female pop for fertility calculation
  fert_idx = pIDX_FERT:(pIDX_FERT + pAG_FERT - 1)
  last_fem_pop = pop[:, FEMALE, fert_idx]

  # Age pop one year
  term_pop = pop[:, :, end]
  pop[:, :, 2:end] = pop[:, :, 1:(end - 1)]
  pop[:, :, end] = pop[:, :, end] + term_pop


  # Age the HIV positive population
  last_hivpop = hivpop
  last_artpop = artpop
  @. hivpop[:, 2:end, :] =  (1-hiv_ag_prob[:, 2:end]) * last_hivpop[:, 2:end, : ] + hiv_ag_prob[:, 1:end-1] * last_hivpop[:, 1:end-1, :]
  if t > par[:tARTstart]
    @. artpop[:, 2:end, :, :] =  (1-hiv_ag_prob[:, 2:end]) * last_artpop[:, 2:end, :, :] + hiv_ag_prob[:, 1:end-1] * last_artpop[:, 1:end-1, :, :]
  end

  # Lagged birth to youngest age group
  entrant_prev = par[:entrantprev]'[t, :]
  # else
  #   entrant_prev = pregprevlag[t-1] * verttrans_lag[t-1] * paedsurv_lag[t-1]
  # end

  if par[:popadjust]
    @. pop[HIVN, :, 1] =  par[:entrantpop]'[t - 1, :] * (1.0-entrant_prev)
    paedsurv = par[:entrantpop]'[t - 1, :] .* entrant_prev
  else
    @. pop[HIVN, :, 1] = par[:birthslag]'[t-1, g] * cumsurv[t-1, g] * (1.0-entrant_prev / paedsurv_lag[t-1]) + cumnetmigr[t-1, g] * (1.0-pregprevlag[t-1] * netmig_hivprob)
    @. paedsurv = par[:birthslag]'[t-1, :] * cumsurv[t-1, g] * entrant_prev + cumnetmigr[t-1, g] * entrant_prev    
  end

  pop[HIVP, :, 1] = paedsurv
  entrantprev_out[t] = sum(pop[HIVP, :, 1]) / sum(pop[:, :, 1])

  first_remain = (1 .- hiv_ag_prob[:, 1]) .* last_hivpop[:, 1, :]
  off_art_entrants = permutedims(par[:paedsurv_cd4dist], (3, 2, 1))[t, :, :] .* ((1.0 .- par[:entrantartcov]'[t, :]) .* paedsurv)
  hivpop[:, 1, :] = first_remain + off_art_entrants
  if t > par[:tARTstart]
    # Remaining in the first age group
    @. artpop[:, 1, :, :] = (1-hiv_ag_prob[:, 1]) * last_artpop[:, 1, :, :]
    # Ageing in on treatment
    artpop[:, 1, :, :] .+= paedsurv .* permutedims(par[:paedsurv_artcd4dist], (4, 3, 2, 1))[t, :, :, :] .* par[:entrantartcov]'[t, :]
  end

  # non-HIV mortality and migration
  hivpop_ha = pop[HIVP, :, :] * age_binner
  qx = 1 .- par[:Sx][t, :, :]
  ndeaths = pop[HIVN, :, :] .* qx
  pop[HIVN, :, :] .-= ndeaths
  hdeaths = pop[HIVP, :, :] .* qx
  pop[HIVP, :, :] .-= hdeaths

  # net migration
  migrate = @. par[:netmigr][t, :, :] * (1+par[:Sx][t, :, :])/2.0 / (pop[HIVN, :, :] + pop[HIVP, :, :])
  pop[HIVN, :, :] .*= 1 .+ migrate
  hmig = migrate .* pop[HIVP, :, :]
  deathsmig = hmig * age_binner
  pop[HIVP, :, :] .+= hmig
  deathmigrate = deathsmig ./ hivpop_ha
  replace!(deathmigrate, NaN=>0.)
  hivpop[:, :, :] .*= 1 .+ deathmigrate
  if(t > par[:tARTstart])
    artpop[:, :, :, :] .*= 1 .+ deathmigrate
  end
  natdeaths[t, :, :] .= ndeaths .+ hdeaths

  # Fertility
  fertile_pop = sum((last_fem_pop + pop[:, FEMALE, fert_idx])./2 .* par[:asfr][:, t]', dims = 1)
  births_by_ha = dropdims(fertile_pop * age_binner[fert_idx, (hIDX_FERT + 1):(hIDX_FERT+hAG_FERT)], dims = 1)
  births = sum(births_by_ha)
  if t + AGE_START < PROJ_YEARS 
    par[:birthslag][:, (t + AGE_START-1), :] = par[:srb][:, t] * births
  end

  return births, births_by_ha, last_hivpop, last_artpop
end

function preg_women_prev_one_step!(par, t, pregprev, pop_ts, pop, hivpop, artpop, births, last_hivpop, last_artpop)
  # Prevalence among pregnant women
  hivbirths = 0.
  for ha = (hIDX_FERT + 1):(hIDX_FERT + hAG_FERT)
    hivn_ha = 0.
    frr_hivpop_ha = 0.
    for a = (hAG_START[ha] + 1):(hAG_START[ha]+hAG_SPAN[ha])
      hivn_ha += (pop_ts[t-1, HIVN, FEMALE, a] + pop[HIVN, FEMALE, a])/2
    end
    for hm = 1:hDS
      frr_hivpop_ha += par[:frr_cd4][t, ha-hIDX_FERT, hm] * (last_hivpop[FEMALE, ha, hm]+hivpop[FEMALE, ha, hm])/2
      if t == par[:tARTstart] 
        for hu = 1:hTS
          frr_hivpop_ha += par[:frr_art][t, ha-hIDX_FERT, hm, hu] * artpop[FEMALE, ha, hm, hu]/2
        end
      elseif t > par[:tARTstart]
        for hu = 1:hTS
          frr_hivpop_ha += par[:frr_art][t, ha-hIDX_FERT, hm, hu] * (last_artpop[FEMALE, ha, hm, hu]+artpop[FEMALE, ha, hm, hu])/2
        end
      end
    end
  end

  pregprev[t] = hivbirths/births
  if (t + AGE_START) < PROJ_YEARS
    pregprevlag[t + AGE_START] = pregprev[t]
  end
end

function adjust_pop!(par, t, pop, hivpop, artpop)
  for g = 1:NG
    a = 1
    for ha = 1:hAG
      popadj_ha = 0.
      hivpop_ha = 0.
      for i = 1:hAG_SPAN[ha]
        hivpop_ha += pop[HIVP, g, a]
        popadjrate_a = popadjust[t, g, a] = par[:targetpop][t, g, a] / (pop[HIVN, g, a] + pop[HIVP, g, a])
        pop[HIVN, g, a] *= popadjrate_a
        hpopadj_a = (popadjrate_a-1.0) * pop[HIVP, g, a]
        popadj_ha += hpopadj_a
        pop[HIVP, g, a] += hpopadj_a
        a += 1
      end
      popadjrate_ha = hivpop_ha > 0 ? popadj_ha / hivpop_ha : 0.0
      for hm = 1:hDS
        hivpop[g, ha, hm] *= 1+popadjrate_ha
        if t >= par[:tARTstart]
          for hu = 1:hTS
            artpop[g, ha, hm, hu] *= 1+popadjrate_ha
          end
        end
      end
    end
  end
end