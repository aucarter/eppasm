
function disease_treatment!(
  par, pop, hivpop, artpop, grad, t, hts, births_by_ha, hivdeaths_ha, everARTelig_idx,
  anyelig_idx, cd4elig_idx, age_binner
)
  # On-ART deaths
  deaths = permutedims(par[:art_mort], (4, 3, 2, 1)) .* artpop
  hivdeaths_ha .+= DT .* dropdims(sum(deaths, dims = (3, 4)), dims = (3, 4))
  aidsdeaths_art[t, :, :, :, :] .+= DT .* deaths
  gradART = -deaths
  # ART time-since initiation progression
  art_prog = artpop[:, :, :, 1:end - 1] ./ h_art_stage_dur
  gradART[:, :, :, 1:end - 1] .-= art_prog
  gradART[:, :, :, 2:end] .+= art_prog
  # ART dropout
  grad .+= par[:art_dropout][t] .* dropdims(sum(artpop, dims = 4), dims = 4)
  gradART .-= par[:art_dropout][t] .* artpop

  Xartelig_15plus, Xart_15plus, expect_mort_artelig15plus, artelig = gettxelig(t, hivpop, par, artpop, gradART, age_binner, pop, births_by_ha, cd4elig_idx)

  ## ART initation
  initart!(hivpop, hts, par, t, anyelig_idx, grad, artinit, gradART, Xartelig_15plus, Xart_15plus, expect_mort_artelig15plus, artelig)

  gradART .= 0.
  artpop = artpop + DT .* gradART
end


function gettxelig(t, hivpop, par, artpop, gradART, age_binner, pop, births_by_ha, cd4elig_idx)
  # Determime number on ART and eligibile
  prop_elig = [(hm >= cd4elig_idx) ? 1.0 : (hm >= hIDX_CD4_350) ? 1.0 - (1.0-specpop_percelig[t])*(1.0-who34percelig) : specpop_percelig[t] for hm = 1:hDS]
  artelig = permutedims(prop_elig .* permutedims(hivpop, (3,1,2)), (2, 3, 1))
  Xartelig_15plus = sum(artelig)
  expect_mort_artelig15plus = sum(par[:cd4_mort] .* artelig)
  Xart_15plus = sum(artpop + DT * gradART)
  # Eligible pregnant women
  if pw_artelig[t] > 0
    p_fert_idx = pIDX_FERT:(pIDX_FERT + pAG_FERT - 1)
    h_fert_idx = (hIDX_FERT + 1):(hIDX_FERT + hAG_FERT)
    neg_women =  (age_binner[p_fert_idx, :]' * pop[HIVN, FEMALE, p_fert_idx])[h_fert_idx]
    pos_no_tx_women = dropdims(sum(par[:frr_cd4][t, :, :] .* hivpop[FEMALE, h_fert_idx, :], dims = 2), dims = 2)
    pos_tx_women = dropdims(sum(par[:frr_art][t, :, :, :] .* artpop[FEMALE, h_fert_idx, :, :], dims = (2, 3)), dims = (2, 3))
    frr_pop_ha =  neg_women .+ pos_no_tx_women .+ pos_tx_women
    pw_elig = births_by_ha .* par[:frr_cd4][t, :, :] .* hivpop[FEMALE, h_fert_idx, :] ./ frr_pop_ha
    artelig[FEMALE, h_fert_idx, :] .+= pw_elig
    Xartelig_15plus += sum(pw_elig)
    expect_mort_artelig15plus += sum(par[:cd4_mort][FEMALE, h_fert_idx, :] .* pw_elig)
  end

  return Xartelig_15plus, Xart_15plus, expect_mort_artelig15plus, artelig
end

function artnuminit(hts, same, perc, art_num, art_num_prev, Xart_15plus, Xartelig_15plus)
  if same 
    artnum_hts = (0.5-DT*(hts+1))*art_num_prev + (DT*(hts+1)+0.5)*art_num
  else
    curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus)
    artnum_hts = curr_coverage + (art_num - curr_coverage) * DT / (0.5-DT*hts)
  end
  if perc
    artnum_hts = artnum_hts * (Xart_15plus + Xartelig_15plus)
  end
  # Only initiate if number exceeds current on ART
  artnum_init = artnum_hts > Xart_15plus ? artnum_hts - Xart_15plus : 0.
  return artnum_init
end

function splitart(artnum_init, artelig, mxweight, Xartelig_15plus, cd4_mort, expect_mort_artelig15plus, hivpop, grad)
  artinit_hahm = artnum_init * artelig * ((1.0 - mxweight)/Xartelig_15plus + mxweight * cd4_mort / expect_mort_artelig15plus)
  if isnan(artinit_hahm)
    artinit_hahm = 0.
  end
  if artinit_hahm > artelig
    artinit_hahm = artelig
  end
  if artinit_hahm > (hivpop + DT * grad)
    artinit_hahm = hivpop + DT * grad
  end
  return artinit_hahm
end

function initart!(hivpop, hts, par, t, anyelig_idx, grad, artinit, gradART, Xartelig_15plus, Xart_15plus, expect_mort_artelig15plus, artelig)
  lag = DT * (hts + 1) < 0.5 ? 1 : 0
  same = par[:art15plus_isperc][:, t-1-lag] == par[:art15plus_isperc][:, t-lag]
  artnum_init = artnuminit.(hts, same, par[:art15plus_isperc][:, t-lag], par[:art15plus_num][:, t - lag], 
    par[:art15plus_num][:, t - 1 - lag], Xart_15plus, Xartelig_15plus)
  # CD4 and age at ART initiation
  if Bool(par[:med_cd4init_input][t])
  elseif par[:art_alloc_method] == 4 # lowest CD4 first
  else # use a mixture of eligibility and expected mortality for initiation distribution
    artinit_hahm = splitart.(
      artnum_init, artelig[:, hIDX_15PLUS:hAG, anyelig_idx:hDS], par[:art_alloc_mxweight][1], 
      Xartelig_15plus, par[:cd4_mort][:, hIDX_15PLUS:hAG, anyelig_idx:hDS],
      expect_mort_artelig15plus, hivpop[:, hIDX_15PLUS:hAG, anyelig_idx:hDS], 
      grad[:, hIDX_15PLUS:hAG, anyelig_idx:hDS]
    )
  end
  grad[:, hIDX_15PLUS:hAG, anyelig_idx:hDS] .-= artinit_hahm / DT
  gradART[:, hIDX_15PLUS:hAG, anyelig_idx:hDS, ART0MOS] .+= artinit_hahm / DT
  artinit[t, :, hIDX_15PLUS:hAG, anyelig_idx:hDS] .+= artinit_hahm 
end