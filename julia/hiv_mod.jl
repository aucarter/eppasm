function hiv_mod_all(p::NamedTuple, x::Dict, t::Int, out::Dict, births_by_ha::Array{Float64, 1})
  for hts = 0:(p.HIVSTEPS_PER_YEAR - 1)
    ts::Int = (t - 2) * p.HIVSTEPS_PER_YEAR + hts + 1
    if p.proj_steps[ts] >= p.tsEpidemicStart
      hiv_mod_one_step!(p, x, t, hts, ts, out, births_by_ha)
    end
  end
end
  
function hiv_mod_one_step!(p::NamedTuple, x::Dict, t::Int, hts::Int, ts::Int, out::Dict, births_by_ha::Array{Float64, 1})
  disease_transmission!(p, x, t, hts, ts, out)
  disease_prog_mort!(p, x, t, out)
  if t >= p.tARTstart
    disease_treatment!(p, x, t, hts, births_by_ha, out)
  end
  x[:hivpop] += p.DT .* x[:grad]
end

function disease_transmission!(p::NamedTuple, x::Dict, t::Int, hts::Int, ts::Int, out::Dict)
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
        t, hts, ts, out)
    end
    out[:prev15to49_ts][ts] = prevcurr
    infections_a_dt = p.DT .* infections_ts
    out[:infections][t, :, :] .+= infections_a_dt
    x[:pop][p.HIVN, :, :] .-=  infections_a_dt
    x[:pop][p.HIVP, :, :] .+=  infections_a_dt
    infections_ha = infections_ts * p.age_binner
    out[:incid15to49][t] += p.DT * sum(infections_ha[p.h_idx_15_49]) 
    x[:grad] .= infections_ha .* p.cd4_initdist
  end
end

function calc_infections_eppspectrum!(p::NamedTuple, x::Dict, r_ts, iota, t::Int, hts::Int, ts::Int, out)
  # Adjust HIV negative and positive population for ptial year time step
  Xhivn_g = dropdims(sum(x[:pop][p.HIVN, :, p.idx_15_49], dims = 2), dims = 2)
  Xhivn = sum(Xhivn_g)
  Xhivn -= sum(x[:pop][p.HIVN, :, p.pIDX_15TO49] .* (1.0 - p.DT*hts))
  Xhivn += sum(x[:pop][p.HIVN, :, p.pIDX_15TO49 + p.pAG_15TO49] .* (1.0 - p.DT*hts))

  Xhivp = sum(x[:pop][p.HIVP, :, p.idx_15_49])
  Xhivp -= sum(x[:pop][p.HIVP, :, p.pIDX_15TO49] .* (1.0 - p.DT*hts))
  Xhivp += sum(x[:pop][p.HIVP, :, p.pIDX_15TO49 + p.pAG_15TO49] .* (1.0 - p.DT*hts))

  prop_include = fill(1., p.NG, length(p.h_idx_15_49))
  prop_include[:, 1] .= 1.0 .- x[:pop][p.HIVP, :, p.hAG_START[(p.hIDX_15TO49 + 1)] + 1] ./  dropdims(sum(x[:pop][p.HIVP, :, (p.hAG_START[(p.hIDX_15TO49 + 1)] + 1):(p.hAG_START[(p.hIDX_15TO49 + 1)]+p.hAG_SPAN[(p.hIDX_15TO49 + 1)])], dims = 2), dims = 2) .* (1.0 - p.DT*hts)
  prop_include[:, end] .= x[:pop][p.HIVP, :, p.hAG_START[(p.hIDX_15TO49 + p.hAG_15TO49)] + 1] ./ dropdims(sum(x[:pop][p.HIVP, :, (p.hAG_START[(p.hIDX_15TO49 + p.hAG_15TO49)]):(p.hAG_START[(p.hIDX_15TO49 + p.hAG_15TO49) + 1]+p.hAG_SPAN[(p.hIDX_15TO49 + p.hAG_15TO49) + 1])], dims = 2), dims = 2) .* (1.0 - p.DT*hts)
  replace!(prop_include, NaN=>0.)
  Xart = sum(prop_include .* x[:artpop])

  Xtot = Xhivn + Xhivp
  prevcurr = Xhivp / Xtot
  @assert !isnan(prevcurr)
  out[:incrate15to49_ts][ts] = r_ts * (Xhivp - Xart + p.relinfectART * Xart) / Xtot + iota

  # incidence by sex
  incrate15to49_g = out[:incrate15to49_ts][ts] * [1, p.incrr_sex[t]] * sum(x[:pop][p.HIVN, :, p.idx_15_49]) / (sum(x[:pop][p.HIVN, 1, p.idx_15_49]) + p.incrr_sex[t]*sum(x[:pop][p.HIVN, 2, p.idx_15_49]))

  # annualized infections by age and sex
  Xhivn_incagerr = dropdims(sum(p.incrr_age[t, :, p.idx_15_49] .* x[:pop][p.HIVN, :, p.idx_15_49], dims = 2), dims = 2)
  infections_ts = zeros(p.NG, p.pAG)
  @. infections_ts = x[:pop][p.HIVN, :, :] * incrate15to49_g * p.incrr_age[t, :, :] * Xhivn_g[:] / Xhivn_incagerr[:]
  return infections_ts, prevcurr
end

 
function disease_prog_mort!(p::NamedTuple, x::Dict, t::Int, out::Dict)
  # Disease progression
  prog = p.cd4_prog .* x[:hivpop][:, :, 1:end - 1]
  x[:grad][:, :, 1:end - 1] -= prog
  x[:grad][:, :, 2:end] += prog 
  
  # Deaths
  x[:grad] .-= p.cd4_mort .* x[:hivpop]
  noart_deaths = p.DT .* p.cd4_mort .* x[:hivpop]
  out[:aidsdeaths_noart][t, :, :, :] .+= noart_deaths
  
  art_deaths = p.DT .* p.art_mort .* x[:artpop]
  out[:aidsdeaths_art][t, :, :, :, :] .+= art_deaths
  
  hivdeaths = dropdims(sum(noart_deaths, dims = 3), dims = 3) .+ 
    dropdims(sum(art_deaths, dims = (3, 4)), dims = (3, 4)) 
  hivpop_props =  x[:pop][p.HIVP, :, :] ./ ((x[:pop][p.HIVP, :, :] * p.age_binner) * p.age_binner')
  hivdeaths_p = hivdeaths * p.age_binner' .* hivpop_props
  # @assert sum(hivdeaths) - sum(hivdeaths_p) < 0.001
  x[:pop][p.HIVP, :, :] .-= hivdeaths_p
  out[:hivdeaths][t, :, :] .+= hivdeaths_p
end