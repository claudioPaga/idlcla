function simulate_calib_al_line, n_counts_in_line, detx, detymin, detymax, sctife, pctife

  DEVICE,DECOMPOSED=0.
  tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1

  al_sigma = 30.0
  al_energy = 1486.0
  fe_energy = 5890.0

                                ; Generate uniformely distributed
                                ; dety coords in the segment of column
                                ; [detymin, detymax]
  
  dety = randomu(seed, n_counts_in_line)*(detymax-detymin)+detymin

                                ; Generate intrinsic energy, using a
                                ; standard deviation of 30eV at Al
                                ; K-alpha

   en_start = randomn(seed, n_counts_in_line)*al_sigma + al_energy
  
   sctial = sctife * (al_energy/fe_energy)^(-0.69623)
   pctial = pctife * (al_energy/fe_energy)^(-0.69623)

   ;print, 'CTI (s/p) at Al = ', sctial, pctial
                                ; Calculate expected damaged energy
                                ; after transfers for each photon
   
   damaged_parallel = en_start * (1.-pctial)^dety
   damaged_serial = damaged_parallel * (1.-sctial)^detx

   losses_electrons = (en_start - damaged_serial)/3.65

   losses_electrons_parallel = (en_start - damaged_parallel)/3.65
   losses_electrons_serial = (damaged_parallel - damaged_serial)/3.65
   
                                ; Randomise the process for each photon
   probability_capture_in_single_transfer = mean(losses_electrons/(dety+detx))

   probability_capture_in_single_transfer_parallel = mean(losses_electrons_parallel/dety)
   probability_capture_in_single_transfer_serial = mean(losses_electrons_serial/detx)
   
   sigma = ((dety+detx)*probability_capture_in_single_transfer * (1.0-probability_capture_in_single_transfer))^0.5
   
   losses_electrons_stocastic = (en_start - damaged_serial)/3.65 + randomu(seed, n_counts_in_line, /norm)*sigma
;Check losses are not negative (can happen in a theoretical stocastic process)
   index_negative = where(losses_electrons_stocastic lt 0., n_negative)
   if n_negative gt 0 then losses_electrons_stocastic[index_negative] = 0.0
  
   damaged_stocastic = en_start - losses_electrons_stocastic * 3.65
 
   plothist, en_start, xhist_start, yhist_start, background=1, color=0, xr=[1000, 1600], bin = 5
   plothist, damaged_parallel, /overplot, color=4, bin = 5
   plothist, damaged_serial, /overplot, color=2, bin = 5
   plothist, damaged_stocastic, xhist_damaged_stocastic, yhist_damage_stocastic, bin = 5, /overplot, color=3

                                ; Second method to derive stocastic
                                ; losses, NOTE - Will use this as I
                                ;                believe it's
                                ;                more representative
                                ;                of what is happening
                                ;                and distinguishes
                                ;                serial and parallel CTIs

   losses_electrons_stocastic_binomial = dblarr(n_counts_in_line)
   for count_photons = 0, n_counts_in_line-1 do losses_electrons_stocastic_binomial[count_photons] = randomu(seed, 1, binomial = [dety[count_photons], probability_capture_in_single_transfer_parallel]) + randomu(seed, 1, binomial = [detx, probability_capture_in_single_transfer_serial])
   damaged_stocastic_binomial = en_start - losses_electrons_stocastic_binomial * 3.65

   plothist, en_start, xhist_start, yhist_start, background=1, color=0, xr=[1000, 1600], bin = 5
   plothist, damaged_stocastic, xhist_damaged_stocastic, yhist_damage_stocastic, bin = 5, /overplot, color=3
   plothist, damaged_stocastic_binomial, xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, bin = 5, /overplot, color=0
      
   res_gaussfit_start = gaussfit(xhist_start, yhist_start, coeff_fit_start, nterms=3, sigma = sigma_start)
   ;print, 'Best fit gauss coeff of original line ', coeff_fit_start
   ;print, 'Best fit sigma of original line: ',sigma_start
   ;stop
   ploterror, xhist_start, yhist_start, sqrt(yhist_start), background=1, color=0, psym=1,  ERRCOLOR = 0
   oplot, xhist_start, res_gaussfit_start, color=3  

; Calculate expected damaged energy after transfers
   

   ploterror, xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, sqrt(yhist_damage_stocastic_bin), background=1, color=0, psym=1,  ERRCOLOR = 0
   ;errors = sqrt(yhist_damage_stocastic)
   res_gaussfit = gaussfit(xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, coeff_fit, nterms=3, sigma = sigma)
   ;print, 'Best fit gauss coeff of damaged line, SERIAL + PARALLEL: ', coeff_fit
   ;print, 'Best fit sigma of damaged line, SERIAL + PARALLEL: ',sigma
   ;print, 'Line centroid: ', coeff_fit[1], sigma[1]
   oplot, xhist_damaged_stocastic_bin, res_gaussfit, color=3  

   cti_measured = 1.-(coeff_fit[1]/al_energy)^(1./(detx+(detymin+detymax)*0.5))
   ;print, 'CTI from this observation = ', cti_measured



   ;PLOTS FOR BOL CASES
   
   if pctife lt  1.E-06 then begin
   
   !P.Multi = [0, 1, 2, 0, 1] 
   
   plothist, en_start, xhist_start, yhist_start, background=1, color=0, xr=[550, 1600], bin = 5, xtit='En (eV)', ytit = 'Counts', tit = 'Black = Model Al K-alpha line, Red = Observed', charsize = 0.8, boxplot = 0, xst=1
   plothist, damaged_stocastic_binomial, xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, bin = 5, /overplot, color=3, boxplot = 0
   ploterror, xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, sqrt(yhist_damage_stocastic_bin), background=1, color=0, psym=1,  ERRCOLOR = 0, xtit='En (eV)', ytit = 'Counts', tit = 'Fit of observed Al K-alpha', charsize = 0.8, xr=[550, 1600], xst=1
   stringx = 'X = '+strtrim(string(detx), 2)
   stringy = 'Y = ['+strtrim(string(detymin), 2)+','+strtrim(string(detymax), 2)+']'
   stringfit = 'En = '+strtrim(string(round(coeff_fit[1])),2)+'+/-'+strtrim(string(round(sigma[1])),2)
   xyouts, 0.15, 0.4, stringx, /norm, color=0
   xyouts, 0.15, 0.37, stringy, /norm, color=0
   xyouts, 0.15, 0.34, stringfit, /norm, color=0
   oplot, xhist_damaged_stocastic_bin, res_gaussfit, color=4     

      
                                ; Postscript plot
   !P.Multi = [0, 1, 2, 0, 1] 
   entry_device = !d.name   
   plotfilenameps = 'calib_al_line_simulation_detx'+strtrim(string(detx),2)+'_detymin'+strtrim(string(detymin),2)+'_detymax'+ strtrim(string(detymax),2)+'_bol.ps'
     set_plot,'ps'
   device, filename = plotfilenameps, /color
   plothist, en_start, xhist_start, yhist_start, background=1, color=0, xr=[550, 1600], bin = 5, xtit='En (eV)', ytit = 'Counts', tit = 'Black = Model Al K-alpha line, Red = Observed', charsize = 0.8, boxplot = 0, xst=1
   plothist, damaged_stocastic_binomial, xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, bin = 5, /overplot, color=3, boxplot = 0
   ploterror, xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, sqrt(yhist_damage_stocastic_bin), background=1, color=0, psym=1,  ERRCOLOR = 0, xtit='En (eV)', ytit = 'Counts', tit = 'Fit of observed Al K-alpha', charsize = 0.8, xr=[550, 1600], xst=1
   stringx = 'X = '+strtrim(string(detx), 2)
   stringy = 'Y = ['+strtrim(string(detymin), 2)+','+strtrim(string(detymax), 2)+']'
   stringfit = 'En (eV)= '+strtrim(string(round(coeff_fit[1])),2)+'+/-'+strtrim(string(round(sigma[1])),2)
   xyouts, 0.12, 0.4, stringx, /norm, color=0, charsize = 0.9
   xyouts, 0.12, 0.37, stringy, /norm, color=0, charsize = 0.9
   xyouts, 0.12, 0.34, stringfit, /norm, color=0, charsize = 0.9
   oplot, xhist_damaged_stocastic_bin, res_gaussfit, color=4     
   device,/close  
   set_plot,entry_device
   !P.Multi = [0, 1, 1, 0, 1] 


endif else begin

   ;PLOTS FOR EOL CASES



!P.Multi = [0, 1, 2, 0, 1] 
   
   plothist, en_start, xhist_start, yhist_start, background=1, color=0, xr=[550, 1600], bin = 5, xtit='En (eV)', ytit = 'Counts', tit = 'Black = Model Al K-alpha line, Red = Observed', charsize = 0.8, boxplot = 0, xst=1
   plothist, damaged_stocastic_binomial, xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, bin = 5, /overplot, color=3, boxplot = 0
   ploterror, xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, sqrt(yhist_damage_stocastic_bin), background=1, color=0, psym=1,  ERRCOLOR = 0, xtit='En (eV)', ytit = 'Counts', tit = 'Fit of observed Al K-alpha', charsize = 0.8, xr=[550, 1600], xst=1
   stringx = 'X = '+strtrim(string(detx), 2)
   stringy = 'Y = ['+strtrim(string(detymin), 2)+','+strtrim(string(detymax), 2)+']'
   stringfit = 'En = '+strtrim(string(round(coeff_fit[1])),2)+'+/-'+strtrim(string(round(sigma[1])),2)
   xyouts, 0.15, 0.4, stringx, /norm, color=0
   xyouts, 0.15, 0.37, stringy, /norm, color=0
   xyouts, 0.15, 0.34, stringfit, /norm, color=0
   oplot, xhist_damaged_stocastic_bin, res_gaussfit, color=4     

      
                                ; Postscript plot
   !P.Multi = [0, 1, 2, 0, 1] 
   entry_device = !d.name   
   plotfilenameps = 'calib_al_line_simulation_detx'+strtrim(string(detx),2)+'_detymin'+strtrim(string(detymin),2)+'_detymax'+ strtrim(string(detymax),2)+'_eol.ps'
     set_plot,'ps'
   device, filename = plotfilenameps, /color
   plothist, en_start, xhist_start, yhist_start, background=1, color=0, xr=[550, 1600], bin = 5, xtit='En (eV)', ytit = 'Counts', tit = 'Black = Model Al K-alpha line, Red = Observed', charsize = 0.8, boxplot = 0
   plothist, damaged_stocastic_binomial, xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, bin = 5, /overplot, color=3, boxplot = 0
   ploterror, xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, sqrt(yhist_damage_stocastic_bin), background=1, color=0, psym=1,  ERRCOLOR = 0, xtit='En (eV)', ytit = 'Counts', tit = 'Fit of observed Al K-alpha', charsize = 0.8, xr=[550, 1600], xst=1
   stringx = 'X = '+strtrim(string(detx), 2)
   stringy = 'Y = ['+strtrim(string(detymin), 2)+','+strtrim(string(detymax), 2)+']'
   stringfit = 'En (eV)= '+strtrim(string(round(coeff_fit[1])),2)+'+/-'+strtrim(string(round(sigma[1])),2)
   xyouts, 0.76, 0.4, stringx, /norm, color=0, charsize = 0.9
   xyouts, 0.76, 0.37, stringy, /norm, color=0, charsize = 0.9
   xyouts, 0.76, 0.34, stringfit, /norm, color=0, charsize = 0.9
   oplot, xhist_damaged_stocastic_bin, res_gaussfit, color=4     
   device,/close  
   set_plot,entry_device
   !P.Multi = [0, 1, 1, 0, 1] 
endelse
array_return = [detx, detymin, detymax, coeff_fit[1], sigma[1], cti_measured]
return, array_return
stop

end
