
function disease_treatment!(
  β, x, t, hts, births_by_ha, hivdeaths_ha, out_dict
)
  cd4elig_idx = β[:artcd4elig_idx][t] 
  anyelig_idx = (β[:specpop_percelig][t] > 0 || β[:pw_artelig][t] > 0) ? 1 : (β[:who34percelig] > 0) ? hIDX_CD4_350 : cd4elig_idx
  β[:everARTelig_idx] = (anyelig_idx < β[:everARTelig_idx]) ? anyelig_idx : β[:everARTelig_idx]
  # On-ART deaths
  deaths = β[:art_mort] .* x[:artpop]
  hivdeaths_ha .+= DT .* dropdims(sum(deaths, dims = (3, 4)), dims = (3, 4))
  out_dict[:aidsdeaths_art][t, :, :, :, :] .+= DT .* deaths
  x[:gradART] = -deaths
  # ART time-since initiation progression
  art_prog = x[:artpop][:, :, :, 1:end - 1] ./ β[:ss][:h_art_stage_dur]
  x[:gradART][:, :, :, 1:end - 1] .-= art_prog
  x[:gradART][:, :, :, 2:end] .+= art_prog
  # ART dropout
  x[:grad] .+= β[:art_dropout][t] .* dropdims(sum(x[:artpop], dims = 4), dims = 4)
  x[:gradART] .-= β[:art_dropout][t] .* x[:artpop]

  # Eligibility
  Xartelig_15plus, Xart_15plus, expect_mort_artelig15plus, artelig = gettxelig(t, x, β, births_by_ha, cd4elig_idx)

  ## ART initation
  initart!(x, hts, β, t, anyelig_idx, Xartelig_15plus, Xart_15plus, expect_mort_artelig15plus, artelig, out_dict)
  x[:artpop] .+= DT .* x[:gradART]
end

function gettxelig(t, x, β,  births_by_ha, cd4elig_idx)
  # Determime number on ART and eligibile
  prop_elig = [(hm >= cd4elig_idx) ? 1.0 : (hm >= hIDX_CD4_350) ? 1.0 - (1.0-β[:specpop_percelig][t])*(1.0-β[:who34percelig]) : β[:specpop_percelig][t] for hm = 1:hDS]
  artelig = permutedims(prop_elig .* permutedims(x[:hivpop], (3,1,2)), (2, 3, 1))
  Xartelig_15plus = dropdims(sum(artelig, dims = (2, 3)), dims = (2, 3))
  expect_mort_artelig15plus = dropdims(sum(β[:cd4_mort] .* artelig, dims = (2, 3)), dims = (2, 3))
  Xart_15plus = dropdims(sum(x[:artpop] + DT * x[:gradART], dims = (2, 3, 4)), dims = (2, 3, 4))
  # Eligible pregnant women
  if β[:pw_artelig][t] > 0
    p_fert_idx = pIDX_FERT:(pIDX_FERT + pAG_FERT - 1)
    h_fert_idx = (hIDX_FERT + 1):(hIDX_FERT + hAG_FERT)
    neg_women =  (β[:age_binner][p_fert_idx, :]' * x[:pop][HIVN, FEMALE, p_fert_idx])[h_fert_idx]
    pos_no_tx_women = dropdims(sum(β[:frr_cd4][t, :, :] .* x[:hivpop][FEMALE, h_fert_idx, :], dims = 2), dims = 2)
    pos_tx_women = dropdims(sum(β[:frr_art][t, :, :, :] .* x[:artpop][FEMALE, h_fert_idx, :, :], dims = (2, 3)), dims = (2, 3))
    frr_pop_ha =  neg_women .+ pos_no_tx_women .+ pos_tx_women
    pw_elig = births_by_ha .* β[:frr_cd4][t, :, :] .* x[:hivpop][FEMALE, h_fert_idx, :] ./ frr_pop_ha
    artelig[FEMALE, h_fert_idx, :] .+= pw_elig
    Xartelig_15plus[FEMALE] += sum(pw_elig)
    expect_mort_artelig15plus[FEMALE] += sum(β[:cd4_mort][FEMALE, h_fert_idx, :] .* pw_elig)
  end

  return Xartelig_15plus, Xart_15plus, expect_mort_artelig15plus, artelig
end


function initart!(x, hts, β, t, anyelig_idx, Xartelig_15plus, Xart_15plus, expect_mort_artelig15plus, artelig, out_dict)
  lag = DT * (hts + 1) < 0.5 ? 1 : 0
  same = β[:art15plus_isperc][:, t-1-lag] == β[:art15plus_isperc][:, t-lag]
  artnum_init = artnuminit.(hts, same, β[:art15plus_isperc][:, t-lag], β[:art15plus_num][t - lag, :], 
    β[:art15plus_num][t - 1 - lag, :], Xart_15plus, Xartelig_15plus, lag)
  # CD4 and age at ART initiation
  if Bool(β[:med_cd4init_input][t])
  elseif β[:art_alloc_method] == 4 # lowest CD4 first
  else # use a mixture of eligibility and expected mortality for initiation distribution
    artinit_weight = @. artelig * ((1.0 - β[:art_alloc_mxweight])/Xartelig_15plus + β[:art_alloc_mxweight] * β[:cd4_mort] / expect_mort_artelig15plus)
    artinit_hahm = artnum_init .* artinit_weight
    artinit_hahm[artinit_hahm .> artelig] .= artelig[artinit_hahm .> artelig]
    expected_hivpop = x[:hivpop] .+ DT .* x[:grad]
    artinit_hahm[artinit_hahm .> expected_hivpop] .= expected_hivpop[artinit_hahm .> expected_hivpop]
  end
  x[:grad] .-= artinit_hahm / DT
  x[:gradART] .+= artinit_hahm / DT
  out_dict[:artinit][t, :, :, :] .+= artinit_hahm 
end

function artnuminit(hts, same, perc, art_num, art_num_prev, Xart_15plus, Xartelig_15plus, lag)
  if same 
    artnum_hts = (1.5 - lag - DT*(hts+1))*art_num_prev + (DT*(hts+1) - 0.5 + lag)*art_num
  else
    curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus)
    artnum_hts = curr_coverage + (art_num - curr_coverage) * DT / (1.5 - lag - DT*hts)
  end
  if perc
    artnum_hts = artnum_hts * (Xart_15plus + Xartelig_15plus)
  end
  # Only initiate if number exceeds current on ART
  artnum_init = artnum_hts > Xart_15plus ? artnum_hts - Xart_15plus : 0.
  return artnum_init
end