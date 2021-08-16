function pop_project_one_step!(β, t, x, age_binner, max_age_binner, out_dict)
  ## Population projection
  # Save last year fertile female pop for fertility calculation
  last_fert_pop = copy(x[:pop][:, FEMALE, fert_idx])

  # Age pop one year
  term_pop = x[:pop][:, :, end]
  x[:pop][:, :, 2:end] = x[:pop][:, :, 1:(end - 1)]
  x[:pop][:, :, end] = x[:pop][:, :, end] + term_pop

  # Calculate probability of aging from one HIV positive age bin to the next
  if sum(x[:pop][HIVP, :, :]) > 0
    hiv_ag_prob = zeros()
    hiv_ag_prob = (x[:pop][HIVP, :, :] * max_age_binner) ./ (x[:pop][HIVP, :, :] * age_binner)
    replace!(hiv_ag_prob, NaN=>0.)
    hiv_ag_prob[:, end] .= 0
  else
    hiv_ag_prob = zeros(NG, hAG)
  end

  # Age the HIV positive population
  last_hivpop = copy(x[:hivpop])
  last_artpop = copy(x[:artpop])
  @. x[:hivpop][:, 2:end, :] =  (1-hiv_ag_prob[:, 2:end]) * last_hivpop[:, 2:end, : ] + hiv_ag_prob[:, 1:end-1] * last_hivpop[:, 1:end-1, :]
  if t > β[:tARTstart]
    @. x[:artpop][:, 2:end, :, :] =  (1-hiv_ag_prob[:, 2:end]) * last_artpop[:, 2:end, :, :] + hiv_ag_prob[:, 1:end-1] * last_artpop[:, 1:end-1, :, :]
  end

  # Lagged birth to youngest age group
  entrant_prev = β[:entrantprev]'[t, :]
  # else
  #   entrant_prev = out_dict[:pregprevlag][t-1] * β[:verttrans_lag][t-1] * β[:paedsurv_lag][t-1]
  # end

  if β[:popadjust]
    @. x[:pop][HIVN, :, 1] =  β[:entrantpop]'[t - 1, :] * (1.0-entrant_prev)
    paedsurv = β[:entrantpop]'[t - 1, :] .* entrant_prev
  else
    @. x[:pop][HIVN, :, 1] = β[:birthslag]'[t-1, g] * cumsurv[t-1, g] * (1.0-entrant_prev / β[:paedsurv_lag][t-1]) + cumnetmigr[t-1, g] * (1.0-out_dict[:pregprevlag][t-1] * netmig_hivprob)
    @. paedsurv = β[:birthslag]'[t-1, :] * cumsurv[t-1, g] * entrant_prev + cumnetmigr[t-1, g] * entrant_prev    
  end

  x[:pop][HIVP, :, 1] = paedsurv
  out_dict[:entrantprev_out][t] = sum(x[:pop][HIVP, :, 1]) / sum(x[:pop][:, :, 1])

  first_remain = (1 .- hiv_ag_prob[:, 1]) .* last_hivpop[:, 1, :]
  off_art_entrants = β[:paedsurv_cd4dist][t, :, :] .* ((1.0 .- β[:entrantartcov]'[t, :]) .* paedsurv)
  x[:hivpop][:, 1, :] = first_remain + off_art_entrants
  if t > β[:tARTstart]
    # Remaining in the first age group
    @. x[:artpop][:, 1, :, :] = (1-hiv_ag_prob[:, 1]) * last_artpop[:, 1, :, :]
    # Ageing in on treatment
    x[:artpop][:, 1, :, :] .+= paedsurv .* β[:paedsurv_artcd4dist][t, :, :, :] .* β[:entrantartcov]'[t, :]
  end

  # non-HIV mortality and migration
  hivpop_ha = x[:pop][HIVP, :, :] * age_binner
  qx = 1 .- β[:Sx][t, :, :]
  ndeaths = x[:pop][HIVN, :, :] .* qx
  x[:pop][HIVN, :, :] .-= ndeaths
  hdeaths = x[:pop][HIVP, :, :] .* qx
  x[:pop][HIVP, :, :] .-= hdeaths
  out_dict[:natdeaths][t, :, :] .= ndeaths .+ hdeaths

  # net migration
  migrate = @. β[:netmigr][t, :, :] * (1+β[:Sx][t, :, :])/2.0 / (x[:pop][HIVN, :, :] + x[:pop][HIVP, :, :])
  x[:pop][HIVN, :, :] .*= 1 .+ migrate
  hmig = migrate .* x[:pop][HIVP, :, :]
  deathsmig = hmig * age_binner
  x[:pop][HIVP, :, :] .+= hmig
  deathmigrate = deathsmig ./ hivpop_ha
  replace!(deathmigrate, NaN=>0.)
  x[:hivpop][:, :, :] .*= 1 .+ deathmigrate
  if(t > β[:tARTstart])
    x[:artpop][:, :, :, :] .*= 1 .+ deathmigrate
  end
  

  # Fertility
  fertile_pop = sum((last_fert_pop + x[:pop][:, FEMALE, fert_idx])./2 .* β[:asfr][:, t]', dims = 1)
  births_by_ha = dropdims(fertile_pop * age_binner[fert_idx, (hIDX_FERT + 1):(hIDX_FERT+hAG_FERT)], dims = 1)
  births = sum(births_by_ha)
  if t + AGE_START < PROJ_YEARS 
    β[:birthslag][:, (t + AGE_START-1), :] = β[:srb][:, t] * births
  end

  return births, births_by_ha, last_hivpop, last_artpop
end

function preg_women_prev_one_step!(β, t, x, births, last_hivpop, last_artpop, out_dict)
  # Prevalence among pregnant women
  hivbirths = 0.
  for ha = (hIDX_FERT + 1):(hIDX_FERT + hAG_FERT)
    hivn_ha = 0.
    frr_hivpop_ha = 0.
    for a = (hAG_START[ha] + 1):(hAG_START[ha]+hAG_SPAN[ha])
      hivn_ha += (out_dict[:pop][t-1, HIVN, FEMALE, a] + x[:pop][HIVN, FEMALE, a])/2
    end
    for hm = 1:hDS
      frr_hivpop_ha += β[:frr_cd4][t, ha-hIDX_FERT, hm] * (last_hivpop[FEMALE, ha, hm]+x[:hivpop][FEMALE, ha, hm])/2
      if t == β[:tARTstart] 
        for hu = 1:hTS
          frr_hivpop_ha += β[:frr_art][t, ha-hIDX_FERT, hm, hu] * x[:artpop][FEMALE, ha, hm, hu]/2
        end
      elseif t > β[:tARTstart]
        for hu = 1:hTS
          frr_hivpop_ha += β[:frr_art][t, ha-hIDX_FERT, hm, hu] * (last_artpop[FEMALE, ha, hm, hu]+x[:artpop][FEMALE, ha, hm, hu])/2
        end
      end
    end
  end

  out_dict[:pregprev][t] = hivbirths/births
  if (t + AGE_START) < PROJ_YEARS
    out_dict[:pregprevlag][t + AGE_START] = hivbirths/births
  end
end

function adjust_pop!(β, t, x, out_dict, age_binner)
  out_dict[:popadjust][t, :, :] .=  β[:targetpop][t, :, :] ./ dropdims(sum(x[:pop], dims = 1), dims = 1)
  x[:pop][HIVN, :, :] .*= out_dict[:popadjust][t, :, :]
  popadj_a = (out_dict[:popadjust][t, :, :] .- 1.0) .* x[:pop][HIVP, :, :]
  popadjrate_ha = (popadj_a * age_binner)./ (x[:pop][HIVP, :, :] * age_binner)
  popadjrate_ha[popadjrate_ha .<= 0] .= 0.
  x[:hivpop] .*= 1 .+ popadjrate_ha
  if t >= β[:tARTstart]
    x[:artpop] .*= 1 .+ popadjrate_ha
  end
end