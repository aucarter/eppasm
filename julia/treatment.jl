
function disease_treatment!(p::NamedTuple, x::Dict, t::Int, hts::Int, births_by_ha::Array{Float64, 1}, out::Dict)
  cd4elig_idx::Int = p.artcd4elig_idx[t] 
  # anyelig_idx = (p.specpop_percelig[t] > 0 || p.pw_artelig[t] > 0) ? 1 : (p.who34percelig > 0) ? p.hIDX_CD4_350 : cd4elig_idx
  # p.everARTelig_idx = (anyelig_idx < p.everARTelig_idx) ? anyelig_idx : p.everARTelig_idx
  # On-ART deaths
  x[:gradART] .= -p.art_mort .* x[:artpop]
  # ART time-since initiation progression
  art_prog = x[:artpop][:, :, :, 1:end - 1] ./ p.ss[:h_art_stage_dur]
  x[:gradART][:, :, :, 1:end - 1] .-= art_prog
  x[:gradART][:, :, :, 2:end] .+= art_prog
  # ART dropout
  x[:grad] .+= p.art_dropout[t] .* sum(x[:artpop], dims = 4)[:, :, :, 1]
  x[:gradART] .-= p.art_dropout[t] .* x[:artpop]

  # Eligibility
  Xartelig_15plus, Xart_15plus, artelig = gettxelig(t, x, p, births_by_ha, cd4elig_idx)

  ## ART initation
  initart!(x, hts, p, t, Xartelig_15plus, Xart_15plus, artelig, out)
  x[:artpop] .+= p.DT .* x[:gradART]
end

function gettxelig(t::Int, x::Dict, p::NamedTuple,  births_by_ha::Array{Float64, 1}, cd4elig_idx::Int)::Tuple{Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 3}}
  # Determime number on ART and eligibile
  prop_elig::Array{Float64, 1} = [(hm >= cd4elig_idx) ? 1.0 : (hm >= p.hIDX_CD4_350) ? 1.0 - (1.0-p.specpop_percelig[t])*(1.0-p.who34percelig) : p.specpop_percelig[t] for hm = 1:p.hDS]
  artelig::Array{Float64, 3} = permutedims(prop_elig .* permutedims(x[:hivpop], (3,1,2)), (2, 3, 1))
  Xartelig_15plus::Array{Float64, 1} = dropdims(sum(artelig, dims = (2, 3)), dims = (2, 3))
  Xart_15plus::Array{Float64, 1} = sum(sum(sum(x[:artpop] + p.DT * x[:gradART], dims = 4), dims = 3), dims = 2)[:, 1, 1, 1]
  # Eligible pregnant women
  if p.pw_artelig[t] > 0
    neg_women =  (p.age_binner[p.fert_idx, :]' * x[:pop][p.HIVN, p.FEMALE, p.fert_idx])[p.h_fert_idx]
    pos_no_tx_women = dropdims(sum(p.frr_cd4[t, :, :] .* x[:hivpop][p.FEMALE, p.h_fert_idx, :], dims = 2), dims = 2)
    pos_tx_women = dropdims(sum(p.frr_art[t, :, :, :] .* x[:artpop][p.FEMALE, p.h_fert_idx, :, :], dims = (2, 3)), dims = (2, 3))
    frr_pop_ha =  neg_women .+ pos_no_tx_women .+ pos_tx_women
    pw_elig = births_by_ha .* p.frr_cd4[t, :, :] .* x[:hivpop][p.FEMALE, p.h_fert_idx, :] ./ frr_pop_ha
    artelig[p.FEMALE, p.h_fert_idx, :] .+= pw_elig
    Xartelig_15plus[p.FEMALE] += sum(pw_elig)
  end

  return Xartelig_15plus, Xart_15plus, artelig
end

function initart!(x::Dict, hts::Int, p::NamedTuple, t::Int, Xartelig_15plus::Array{Float64, 1}, Xart_15plus::Array{Float64, 1}, artelig::Array{Float64, 3}, out::Dict)
  lag = p.DT * (hts + 1) < 0.5 ? 1 : 0
  same = p.art15plus_isperc[:, t-1-lag] == p.art15plus_isperc[:, t-lag]
  artnum_init = artnuminit.(hts, same, p.art15plus_isperc[:, t-lag], p.art15plus_num[t - lag, :], 
    p.art15plus_num[t - 1 - lag, :], Xart_15plus, Xartelig_15plus, lag, p.DT)
  # CD4 and age at ART initiation
  if Bool(p.med_cd4init_input[t])
  elseif p.art_alloc_method == 4 # lowest CD4 first
  else # use a mixture of eligibility and expected mortality for initiation CD4 distribution
    expect_mort_artelig15plus = dropdims(sum(p.cd4_mort .* artelig, dims = (2, 3)), dims = (2, 3))
    expect_mort_weight = p.cd4_mort ./ expect_mort_artelig15plus
    artinit_weight = @. ((1.0 - p.art_alloc_mxweight)/Xartelig_15plus + p.art_alloc_mxweight * expect_mort_weight)
    # global flag
    # if flag
    #   display(x[:hivpop][1, :, :])
    #   flag = false
    # end
    artinit_hahm = artnum_init .* (artinit_weight .* artelig)
    artinit_hahm[artinit_hahm .> artelig] .= artelig[artinit_hahm .> artelig]
    expected_hivpop = x[:hivpop] .+ p.DT .* x[:grad]
    artinit_hahm[artinit_hahm .> expected_hivpop] .= expected_hivpop[artinit_hahm .> expected_hivpop]
  end
  x[:grad] .-= artinit_hahm / p.DT
  x[:gradART] .+= artinit_hahm / p.DT
  out[:artinit][t, :, :, :] .+= artinit_hahm 
end

function artnuminit(hts::Int, same::Bool, perc::Bool, art_num::Float64, art_num_prev::Float64, Xart_15plus::Float64, Xartelig_15plus::Float64, lag::Int, DT::Float64)::Float64
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