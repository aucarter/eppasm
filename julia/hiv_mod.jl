  # Disease progression and mortality
function disease_prog_mort!(hivdeaths_ha, grad, aidsdeaths_noart, t, everARTelig_idx)
  for g = 1:NG
    for ha = 1:hAG
      for hm = 1:hDS
        cd4mx_scale = 1.
        if scale_cd4_mort & t >= t_ART_start & hm >= everARTelig_idx
          artpop_hahm = 0.
          for hu = 1:hTS
            artpop_hahm += artpop[t, g, ha , hm, hu]
            cd4mx_scale = hivpop[t, g, ha, hm] / (hivpop[t, g, ha, hm] + artpop_hahm)
            deaths = cd4mx_scale * cd4_mort[g, ha, hm] * hivpop[t, g, ha, hm]
            hivdeaths_ha[g, ha] += DT*deaths
            aidsdeaths_noart[t, g, ha, hm] += DT*deaths
            grad[g, ha, hm] = -deaths
            if hm > 1
              grad[g, ha, hm-1] -= cd4_prog[g, ha, hm-1] * hivpop[t, g, ha, hm-1]
              grad[g, ha, hm] += cd4_prog[g, ha, hm-1] * hivpop[t, g, ha, hm-1]            
            end
          end
        end
      end
    end
  end
end

function hiv_mod_one_step(t, hts, anyelig_idx, cd4elig_idx, artpop, hivpop, aidsdeaths_noart, everARTelig_idx)
  ts = (t - 2) * HIVSTEPS_PER_YEAR + hts + 1
  hivdeaths_ha = zeros(NG, hAG)
  grad = zeros(NG, hAG, hDS)
  disease_prog_mort!(hivdeaths_ha, grad, aidsdeaths_noart, t, everARTelig_idx)
  
  if eppmod != EPP_DIRECTINCID
    if eppmod == EPP_RSPLINE
      rvec[ts] = rspline_rvec[ts]
    else
      # Need to implement calc_rtrend_rt
      nothing
    end

    if incidmod == INCIDMOD_EPPSPEC 
      infections_ts, prevcurr = calc_infections_eppspectrum(
        pop, hivpop, artpop, rvec[ts], relinfectART, (projsteps[ts] == tsEpidemicStart) ? iota : 0.0,
        incrr_sex, incrr_age,  t_ART_start, DT, t, hts, hAG_START, hAG_SPAN, ts)
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
          pop[t, HIVN, g, a] -= infections_a_dt
          pop[t, HIVP, g, a] += infections_a_dt
          a += 1
        end

        if ha < hIDX_15TO49+hAG_15TO49 
          incid15to49[t] += DT*infections_ha
        end
          
        for hm = 1:hDS
          grad[g, ha, hm] += infections_ha * cd4_initdist[g, ha, hm]
        end
      end
    end
  end
  # ART progression, mortality, and initiation
  if t >= t_ART_start
    gradART = zeros(NG, hAG, hDS, hTS)
    for g = 1:NG
      for ha = 1:hAG
        for hm = 1:hDS
          for hu = 1:hTS
            deaths = art_mort[g, ha, hm, hu] * artmx_timerr[t, hu] * artpop[t, g, ha, hm, hu]
            hivdeaths_ha[g, ha] += DT*deaths
            aidsdeaths_art[t, g, ha, hm, hu] += DT*deaths
            gradART[g, ha, hm, hu] = -deaths
            if hu < hTS
              gradART[g, ha, hm, hu] += -artpop[t, g, ha, hm, hu] / h_art_stage_dur[hu]
              gradART[g, ha, hm, hu+1] += artpop[t, g, ha, hm, hu] / h_art_stage_dur[hu]
            end     
            # ART dropout
            if art_dropout[t] > 0  
              grad[g, ha, hm] += art_dropout[t] * artpop[t, g, ha, hm, hu]
              gradART[g, ha, hm, hu] -= art_dropout[t] * artpop[t, g, ha, hm, hu]
            end
          end
        end
      end
    end

    # ART initation
    for g = 1:NG
      artelig_hahm = zeros(hAG_15PLUS, hDS)
      Xart_15plus = 0.0
      Xartelig_15plus = 0.0
      expect_mort_artelig15plus = 0.0
      for ha = hIDX_15PLUS:hAG
        for hm = everARTelig_idx:hDS
          if hm >= anyelig_idx
            prop_elig = (hm >= cd4elig_idx) ? 1.0 : (hm >= hIDX_CD4_350) ? 1.0 - (1.0-specpop_percelig[t])*(1.0-who34percelig) : specpop_percelig[t]
            Xartelig_15plus += artelig_hahm[ha-hIDX_15PLUS + 1, hm] = prop_elig * hivpop[t, g, ha, hm] 
            expect_mort_artelig15plus += cd4_mort[g, ha, hm] * artelig_hahm[ha-hIDX_15PLUS + 1, hm]
          end
          for hu = 1:hTS
            Xart_15plus += artpop[t, g, ha, hm, hu] + DT * gradART[g, ha, hm, hu]
          end
        end
        if g == FEMALE && pw_artelig[t] > 0 && ha < hAG_FERT
          frr_pop_ha = 0.
          for a = (hAG_START[ha] + 1):(hAG_START[ha]+hAG_SPAN[ha])
            frr_pop_ha += pop[t, HIVN, g, a]
          end
          for hm = 1:hDS
            frr_pop_ha += frr_cd4[t, ha-hIDX_FERT, hm] * hivpop[t, g, ha, hm]
            for hu = 1:hTS
              frr_pop_ha += frr_art[t, ha-hIDX_FERT, hm, hu] * artpop[t, g, ha, hm, hu]
            end
          end
          for hm = anyelig_idx:cd4elig_idx
            pw_elig_hahm = births_by_ha[ha-hIDX_FERT] * frr_cd4[t, ha-hIDX_FERT, hm] * hivpop[t, g, ha, hm] / frr_pop_ha
            artelig_hahm[ha-hIDX_15PLUS + 1, hm] += pw_elig_hahm
            Xartelig_15plus += pw_elig_hahm
            expect_mort_artelig15plus += cd4_mort[g, ha, hm] * pw_elig_hahm
          end
        end
      end

      # Calculate number on ART at end of ts, based on number or percent
      artnum_hts = 0.0
      if DT * (hts + 1) < 0.5
        if !art15plus_isperc[t-2, g] & !art15plus_isperc[t-1, g]
          artnum_hts = (0.5-DT*(hts+1))*artnum15plus[t-2, g] + (DT*(hts+1)+0.5)*artnum15plus[t-1, g]
        elseif art15plus_isperc[t-2, g] & art15plus_isperc[t-1, g]
          artcov_hts = (0.5-DT*(hts+1))*artnum15plus[t-2, g] + (DT*(hts+1)+0.5)*artnum15plus[t-1, g]
          artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus)
        elseif !art15plus_isperc[t-2, g] & art15plus_isperc[t-1, g]
          curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus)
          artcov_hts = curr_coverage + (artnum15plus[t-1, g] - curr_coverage) * DT / (0.5-DT*hts)
          artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus)
        end
      else
        if !art15plus_isperc[t-1, g] & !art15plus_isperc[t, g]
          artnum_hts = (0.5-DT*(hts+1))*artnum15plus[t-1, g] + (DT*(hts+1)+0.5)*artnum15plus[t, g]
        elseif art15plus_isperc[t-1, g] & art15plus_isperc[t, g]
          artcov_hts = (0.5-DT*(hts+1))*artnum15plus[t-1, g] + (DT*(hts+1)+0.5)*artnum15plus[t, g]
          artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus)
        elseif !art15plus_isperc[t-1, g] & art15plus_isperc[t, g]
          curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus)
          artcov_hts = curr_coverage + (artnum15plus[t, g] - curr_coverage) * DT / (0.5-DT*hts)
          artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus)
        end         
      end
      artinit_hts = artnum_hts > Xart_15plus ? artnum_hts - Xart_15plus : 0

      # median CD4 at initiation inputs
      if Bool(med_cd4init_input[t])
      elseif art_alloc_method == 4 # lowest CD4 first
      else # use a mixture of eligibility and expected mortality for initiation distribution
        for ha = hIDX_15PLUS:hAG
          for hm = anyelig_idx:hDS
            artinit_hahm = artinit_hts * artelig_hahm[ha-hIDX_15PLUS + 1, hm] * ((1.0 - art_alloc_mxweight)/Xartelig_15plus + art_alloc_mxweight * cd4_mort[g, ha, hm] / expect_mort_artelig15plus)
            if isnan(artinit_hahm)
              artinit_hahm = 0.
            end
            if artinit_hahm > artelig_hahm[ha-hIDX_15PLUS + 1, hm]
              artinit_hahm = artelig_hahm[ha-hIDX_15PLUS + 1, hm]
            end
            if artinit_hahm > (hivpop[t, g, ha, hm] + DT * grad[g, ha, hm])
              artinit_hahm = hivpop[t, g, ha, hm] + DT * grad[g, ha, hm]
            end
            grad[g, ha, hm] -= artinit_hahm / DT
            gradART[g, ha, hm, ART0MOS] += artinit_hahm / DT
            artinit[t, g, ha, hm] += artinit_hahm 
          end
        end
      end
    end

    for g = 1:NG
      for ha = 1:hAG
        for hm = everARTelig_idx:hDS
          for hu = 1:hTS
            if isnan(gradART[g, ha, hm, hu])
              error("NaN gradART")
            end
            artpop[t, g, ha, hm, hu] += DT*gradART[g, ha, hm, hu]
          end
        end
      end
    end
  end
  for g = 1:NG
    for ha = 1:hAG
      for hm = everARTelig_idx:hDS
        if isnan(grad[g, ha, hm])
          error("NaN grad")
        end
        hivpop[t, g, ha, hm] += DT*grad[g, ha, hm]
      end
    end
  end

  for g = 1:NG
    hivpop_ha = zeros(hAG)
    a = 1
    for ha = 1:hAG
      hivpop_ha[ha] = 0.0
      for i = 1:hAG_SPAN[ha]
        hivpop_ha[ha] += pop[t, HIVP, g, a]
        a += 1
      end
    end

    a = 1
    for ha = 1:hAG
      if hivpop_ha[ha] > 0
        hivqx_ha = hivdeaths_ha[g, ha] / hivpop_ha[ha]
        for i = 1:hAG_SPAN[ha]
          hivdeaths[t, g, a] += pop[t, HIVP, g, a] * hivqx_ha
          pop[t, HIVP, g, a] *= (1.0-hivqx_ha)
          a += 1
        end
      else
        a += hAG_SPAN[ha]
      end
    end
  end
end