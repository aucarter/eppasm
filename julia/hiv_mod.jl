  # Disease progression and mortality
function disease_prog_mort!(mod_data, hivpop, hivdeaths_ha, grad, aidsdeaths_noart, t, everARTelig_idx)
  prog = mod_data[:cd4_prog] .* hivpop[:, :, 1:end - 1]
  grad[:, :, 1:end - 1] -= prog
  grad[:, :, 2:end] += prog 
  deaths = mod_data[:cd4_mort] .* hivpop
  aidsdeaths_noart[t, :, :, :] .+= DT .* deaths
  hivdeaths_ha .+= DT .* dropdims(sum(deaths, dims = 3), dims = 3)
  grad = grad - deaths
end
function disease_transmission!(mod_data, pop, hivpop, artpop, grad, t, hts, ts)
  if eppmod != EPP_DIRECTINCID
    if eppmod == EPP_RSPLINE
      nothing
    else
      # Need to implement calc_rtrend_rt
      nothing
    end

    if incidmod == INCIDMOD_EPPSPEC 
      infections_ts, prevcurr = calc_infections_eppspectrum!(
        mod_data, pop, hivpop, artpop, rvec[ts], (mod_data[:proj_steps][ts] == mod_data[:tsEpidemicStart]) ? mod_data[:iota] : 0.0,
        t, hts, ts)
    end
    global prev15to49_ts[ts] = prevcurr

    for g = 1:NG
      a = 1
      for ha = 1:hAG
        infections_ha = 0
        for i = 1:hAG_SPAN[ha]
          infections_a = infections_ts[g, a]
          infections_ha += infections_a
          infections_a_dt = DT*infections_a
          infections[t, g, a] += infections_a_dt
          if t == 30
            nothing
          end
          pop[HIVN, g, a] -= infections_a_dt
          pop[HIVP, g, a] += infections_a_dt
          a += 1
        end

        if ha < hIDX_15TO49+hAG_15TO49 
          incid15to49[t] += DT*infections_ha
        end
          
        for hm = 1:hDS
          grad[g, ha, hm] += infections_ha * mod_data[:cd4_initdist][hm, ha, g]
        end
      end
    end
  end
end

function disease_treatment!(
  mod_data, pop, hivpop, artpop, grad, t, hts, births_by_ha, hivdeaths_ha, everARTelig_idx,
  anyelig_idx, cd4elig_idx
)
  # On-ART deaths
  deaths = permutedims(mod_data[:art_mort], (4, 3, 2, 1)) .* artpop
  hivdeaths_ha .+= DT .* dropdims(sum(deaths, dims = (3, 4)), dims = (3, 4))
  aidsdeaths_art[t, :, :, :, :] .+= DT .* deaths
  gradART = -deaths
  # ART time-since initiation progression
  art_prog = artpop[:, :, :, 1:end - 1] ./ h_art_stage_dur
  gradART[:, :, :, 1:end - 1] .-= art_prog
  gradART[:, :, :, 2:end] .+= art_prog
  # ART dropout
  grad .+= mod_data[:art_dropout][t] .* dropdims(sum(artpop, dims = 4), dims = 4)
  gradART .-= mod_data[:art_dropout][t] .* artpop

  # ART initation
  # TODO: Pick up here
  for g = 1:NG
    artelig_hahm = zeros(hAG_15PLUS, hDS)
    Xart_15plus = 0.0
    Xartelig_15plus = 0.0
    expect_mort_artelig15plus = 0.0
    for ha = hIDX_15PLUS:hAG
      for hm = everARTelig_idx:hDS
        if hm >= anyelig_idx
          prop_elig = (hm >= cd4elig_idx) ? 1.0 : (hm >= hIDX_CD4_350) ? 1.0 - (1.0-specpop_percelig[t])*(1.0-who34percelig) : specpop_percelig[t]
          Xartelig_15plus += artelig_hahm[ha-hIDX_15PLUS + 1, hm] = prop_elig * hivpop[g, ha, hm] 
          expect_mort_artelig15plus += mod_data[:cd4_mort][g, ha, hm] * artelig_hahm[ha-hIDX_15PLUS + 1, hm]
        end
        for hu = 1:hTS
          Xart_15plus += artpop[g, ha, hm, hu] + DT * gradART[g, ha, hm, hu]
        end
      end
      if g == FEMALE && pw_artelig[t] > 0 && ha < hAG_FERT
        frr_pop_ha = 0.
        for a = (hAG_START[ha] + 1):(hAG_START[ha]+hAG_SPAN[ha])
          frr_pop_ha += pop[HIVN, g, a]
        end
        for hm = 1:hDS
          frr_pop_ha += mod_data[:frr_cd4][hm, ha-hIDX_FERT, t] * hivpop[g, ha, hm]
          for hu = 1:hTS
            frr_pop_ha += mod_data[:frr_art][hu, hm, ha-hIDX_FERT, t] * artpop[g, ha, hm, hu]
          end
        end
        for hm = anyelig_idx:cd4elig_idx
          pw_elig_hahm = births_by_ha[ha-hIDX_FERT] * mod_data[:frr_cd4][hm, ha-hIDX_FERT, t] * hivpop[g, ha, hm] / frr_pop_ha
          artelig_hahm[ha-hIDX_15PLUS + 1, hm] += pw_elig_hahm
          Xartelig_15plus += pw_elig_hahm
          expect_mort_artelig15plus += mod_data[:cd4_mort][g, ha, hm] * pw_elig_hahm
        end
      end
    end

    # Calculate number on ART at end of ts, based on number or percent
    artnum_hts = 0.0
    if DT * (hts + 1) < 0.5
      if !mod_data[:art15plus_isperc][g, t-2] & !mod_data[:art15plus_isperc][g, t-1]
        artnum_hts = (0.5-DT*(hts+1))*mod_data[:art15plus_num][g, t - 2] + (DT*(hts+1)+0.5)*mod_data[:art15plus_num][g, t - 1]
      elseif mod_data[:art15plus_isperc][g, t-2] & mod_data[:art15plus_isperc][g, t-1]
        artcov_hts = (0.5-DT*(hts+1))*mod_data[:art15plus_num][g, t - 2] + (DT*(hts+1)+0.5)*mod_data[:art15plus_num][g, t - 1]
        artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus)
      elseif !mod_data[:art15plus_isperc][g, t-2] & mod_data[:art15plus_isperc][g, t-1]
        curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus)
        artcov_hts = curr_coverage + (mod_data[:art15plus_num][g, t - 1] - curr_coverage) * DT / (0.5-DT*hts)
        artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus)
      end
    else
      if !mod_data[:art15plus_isperc][g, t-1] & !mod_data[:art15plus_isperc][g, t]
        artnum_hts = (0.5-DT*(hts+1))*mod_data[:art15plus_num][g, t - 1] + (DT*(hts+1)+0.5)*mod_data[:art15plus_num][g, t]
      elseif mod_data[:art15plus_isperc][g, t-1] & mod_data[:art15plus_isperc][g, t]
        artcov_hts = (0.5-DT*(hts+1))*mod_data[:art15plus_num][g, t - 1] + (DT*(hts+1)+0.5)*mod_data[:art15plus_num][g, t]
        artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus)
      elseif !mod_data[:art15plus_isperc][g, t-1] & mod_data[:art15plus_isperc][g, t]
        curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus)
        artcov_hts = curr_coverage + (mod_data[:art15plus_num][g, t] - curr_coverage) * DT / (0.5-DT*hts)
        artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus)
      end         
    end
    artinit_hts = artnum_hts > Xart_15plus ? artnum_hts - Xart_15plus : 0

    # median CD4 at initiation inputs
    if Bool(mod_data[:med_cd4init_input][t])
    elseif mod_data[:art_alloc_method] == 4 # lowest CD4 first
    else # use a mixture of eligibility and expected mortality for initiation distribution
      for ha = hIDX_15PLUS:hAG
        for hm = anyelig_idx:hDS
          artinit_hahm = artinit_hts * artelig_hahm[ha-hIDX_15PLUS + 1, hm] * ((1.0 - mod_data[:art_alloc_mxweight][1])/Xartelig_15plus + mod_data[:art_alloc_mxweight][1] * mod_data[:cd4_mort][g, ha, hm] / expect_mort_artelig15plus)
          if isnan(artinit_hahm)
            artinit_hahm = 0.
          end
          if artinit_hahm > artelig_hahm[ha-hIDX_15PLUS + 1, hm]
            artinit_hahm = artelig_hahm[ha-hIDX_15PLUS + 1, hm]
          end
          if artinit_hahm > (hivpop[g, ha, hm] + DT * grad[g, ha, hm])
            artinit_hahm = hivpop[g, ha, hm] + DT * grad[g, ha, hm]
          end
          grad[g, ha, hm] -= artinit_hahm / DT
          gradART[g, ha, hm, ART0MOS] += artinit_hahm / DT
          artinit[t, g, ha, hm] += artinit_hahm 
        end
      end
    end
  end

  artpop = artpop + DT .* gradART
end


function hiv_mod_one_step!(mod_data, t, hts, anyelig_idx, cd4elig_idx, artpop, hivpop, aidsdeaths_noart, everARTelig_idx, births_by_ha, pop, age_binner)
  
  ts = (t - 2) * HIVSTEPS_PER_YEAR + hts + 1
  hivdeaths_ha = zeros(NG, hAG)
  grad = zeros(NG, hAG, hDS)
  
  disease_prog_mort!(mod_data, hivpop, hivdeaths_ha, grad, aidsdeaths_noart, t, everARTelig_idx)
  disease_transmission!(mod_data, pop, hivpop, artpop, grad, t, hts, ts)
  if t >= mod_data[:tARTstart]
    disease_treatment!(
      mod_data, pop, hivpop, artpop, grad, t, hts, births_by_ha, hivdeaths_ha, everARTelig_idx,
      anyelig_idx, cd4elig_idx
    )
  end

  hivpop = hivpop + DT .* grad
  hivpop_ha = pop[HIVP, :, :] * age_binner
  hivqx_ha = hivdeaths_ha ./ hivpop_ha
  replace!(hivqx_ha, NaN=>0.)
  hivqx_all_age = (hivqx_ha * age_binner')
  hivdeaths[t, :, :] = pop[HIVP, :, :] .*  hivqx_all_age
  pop[HIVP, :, :] .*= 1.0 .- hivqx_all_age
end