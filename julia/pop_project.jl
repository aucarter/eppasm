function pop_project_one_step!(p::NamedTuple, t::Int, x::Dict, out::Dict)::Array{Float64, 1}
  ## Population projection
  # Save last year fertile female pop for fertility calculation
  last_fert_pop = copy(x[:pop][:, p.FEMALE, p.fert_idx])

  # Calculate probability of aging from one HIV positive age bin to the next
  if sum(x[:pop][p.HIVP, :, :]) > 0
    hiv_ag_prob = zeros()
    hiv_ag_prob = (x[:pop][p.HIVP, :, :] * p.max_age_binner) ./ (x[:pop][p.HIVP, :, :] * p.age_binner)
    replace!(hiv_ag_prob, NaN=>0.)
    hiv_ag_prob[:, end] .= 0
  else
    hiv_ag_prob = zeros(p.NG, p.hAG)
  end

  # Age pop one year
  term_pop = x[:pop][:, :, end]
  x[:pop][:, :, 2:end] = x[:pop][:, :, 1:(end - 1)]
  x[:pop][:, :, end] = x[:pop][:, :, end] + term_pop

  # Age the HIV positive population
  @. x[:hivpop][:, 2:end, :] .+=  hiv_ag_prob[:, 1:end-1] * out[:hivpop][t-1, :, 1:end-1, :]
  @. x[:hivpop][:, 1:end-1, :] .-=  hiv_ag_prob[:, 1:end-1] * out[:hivpop][t-1, :, 1:end-1, :]
  if t > p.tARTstart
    @. x[:artpop][:, 2:end, :, :] .+=  hiv_ag_prob[:, 1:end-1] * out[:artpop][t-1, :, 1:end-1, :, :]
    @. x[:artpop][:, 1:end-1, :, :] .-=  hiv_ag_prob[:, 1:end-1] * out[:artpop][t-1, :, 1:end-1, :, :]
  end

  # Lagged birth to youngest age group
  entrant_prev = p.entrantprev[t, :]
  # else
  #   entrant_prev = out[:pregprevlag][t-1] * p.verttrans_lag[t-1] * p.paedsurv_lag[t-1]
  # end

  if p.popadjust
    @. x[:pop][p.HIVN, :, 1] =  p.entrantpop[t - 1, :] * (1.0-entrant_prev)
    paedsurv = p.entrantpop[t - 1, :] .* entrant_prev
  else
    @. x[:pop][p.HIVN, :, 1] = p.birthslag[t-1, g] * cumsurv[t-1, g] * (1.0-entrant_prev / p.paedsurv_lag[t-1]) + cumnetmigr[t-1, g] * (1.0-out[:pregprevlag][t-1] * netmig_hivprob)
    @. paedsurv = p.birthslag[t-1, :] * cumsurv[t-1, g] * entrant_prev + cumnetmigr[t-1, g] * entrant_prev    
  end

  x[:pop][p.HIVP, :, 1] = paedsurv
  out[:entrantprev_out][t] = sum(x[:pop][p.HIVP, :, 1]) / sum(x[:pop][:, :, 1])

  first_remain = (1 .- hiv_ag_prob[:, 1]) .* out[:hivpop][t-1, :, 1, :]
  off_art_entrants = p.paedsurv_cd4dist[t, :, :] .* ((1.0 .- p.entrantartcov[t, :]) .* paedsurv)
  x[:hivpop][:, 1, :] = first_remain + off_art_entrants
  if t > p.tARTstart
    # Remaining in the first age group
    @. x[:artpop][:, 1, :, :] = (1-hiv_ag_prob[:, 1]) * out[:artpop][t-1, :, 1, :, :]
    # Ageing in on treatment
    x[:artpop][:, 1, :, :] .+= paedsurv .* p.paedsurv_artcd4dist[t, :, :, :] .* p.entrantartcov[t, :]
  end

  # non-HIV mortality and migration
  hivpop_ha = x[:pop][p.HIVP, :, :] * p.age_binner
  qx = 1 .- p.Sx[t, :, :]
  ndeaths = x[:pop][p.HIVN, :, :] .* qx
  x[:pop][p.HIVN, :, :] .-= ndeaths
  hdeaths = x[:pop][p.HIVP, :, :] .* qx
  x[:pop][p.HIVP, :, :] .-= hdeaths
  out[:natdeaths][t, :, :] .= ndeaths .+ hdeaths
  
  # net migration
  migrate = @. p.netmigr[t, :, :] * (1+p.Sx[t, :, :])/2.0 / (x[:pop][p.HIVN, :, :] + x[:pop][p.HIVP, :, :])
  x[:pop][p.HIVN, :, :] .*= 1 .+ migrate
  hmig = migrate .* x[:pop][p.HIVP, :, :]
  deathsmig = hmig * p.age_binner
  x[:pop][p.HIVP, :, :] .+= hmig
  deathmigrate = deathsmig ./ hivpop_ha
  replace!(deathmigrate, NaN=>0.)
  x[:hivpop][:, :, :] .*= 1 .+ deathmigrate
  if(t > p.tARTstart)
    x[:artpop][:, :, :, :] .*= 1 .+ deathmigrate
  end
  

  # Fertility
  fertile_pop = sum((last_fert_pop + x[:pop][:, p.FEMALE, p.fert_idx])./2 .* p.asfr[t, :]', dims = 1)
  births_by_ha = dropdims(fertile_pop * p.age_binner[p.fert_idx, (p.hIDX_FERT + 1):(p.hIDX_FERT+p.hAG_FERT)], dims = 1)
  if t + p.AGE_START < p.PROJ_YEARS 
    p.birthslag[(t + p.AGE_START-1), :] = p.srb[t, :] * sum(births_by_ha)
  end

  return births_by_ha
end

function preg_women_prev_one_step!(p::NamedTuple, t::Int, x::Dict, out::Dict, births_by_ha::Array{Float64, 1})
  # Prevalence among pregnant women
  hivbirths = 0.
  for ha = (p.hIDX_FERT + 1):(p.hIDX_FERT + p.hAG_FERT)
    hivn_ha = 0.
    frr_hivpop_ha = 0.
    for a = (p.hAG_START[ha] + 1):(p.hAG_START[ha]+p.hAG_SPAN[ha])
      hivn_ha += (out[:pop][t-1, p.HIVN, p.FEMALE, a] + x[:pop][p.HIVN, p.FEMALE, a])/2
    end
    for hm = 1:p.hDS
      frr_hivpop_ha += p.frr_cd4[t, ha-p.hIDX_FERT, hm] * (out[:hivpop][t-1, p.FEMALE, ha, hm]+x[:hivpop][p.FEMALE, ha, hm])/2
      if t == p.tARTstart 
        for hu = 1:p.hTS
          frr_hivpop_ha += p.frr_art[t, ha-p.hIDX_FERT, hm, hu] * x[:artpop][p.FEMALE, ha, hm, hu]/2
        end
      elseif t > p.tARTstart
        for hu = 1:p.hTS
          frr_hivpop_ha += p.frr_art[t, ha-p.hIDX_FERT, hm, hu] * (out[:artpop][t-1, p.FEMALE, ha, hm, hu]+x[:artpop][p.FEMALE, ha, hm, hu])/2
        end
      end
    end
  end

  out[:pregprev][t] = hivbirths/sum(births_by_ha)
  if (t + p.AGE_START) < p.PROJ_YEARS
    out[:pregprevlag][t + p.AGE_START] = hivbirths/sum(births_by_ha)
  end
end

function adjust_pop!(p::NamedTuple, t::Int, x::Dict, out::Dict)
  out[:popadjust][t, :, :] .=  p.targetpop[t, :, :] ./ dropdims(sum(x[:pop], dims = 1), dims = 1)
  x[:pop][p.HIVN, :, :] .*= out[:popadjust][t, :, :]
  popadj_a = (out[:popadjust][t, :, :] .- 1.0) .* x[:pop][p.HIVP, :, :]
  popadjrate_ha = (popadj_a * p.age_binner)./ (x[:pop][p.HIVP, :, :] * p.age_binner)
  x[:pop][p.HIVP, :, :] .+= popadj_a
  # popadjrate_ha[popadjrate_ha .<= 0] .= 0.
  replace!(popadjrate_ha, NaN=>0.)
  x[:hivpop] .*= 1 .+ popadjrate_ha
  if t >= p.tARTstart
    x[:artpop] .*= 1 .+ popadjrate_ha
  end
end