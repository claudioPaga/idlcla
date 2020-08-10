function simulate_calib_source, n_counts_in_line, detx, detymin, detymax, sctife, pctife

;;;
;;; CP, 6 June 2019
;;; Summary - Simulate CTI distortion of Al K-alpha line
;;;
;;; KEYWORD PARAMETERS:
;;; n_counts_in_line 
;;; detx = Column X coordinate with Al line events
;;; detymin, detymax = Y min, max coordinates of column segment
;;; sctife, pctife = serial and parallel CTI values at Fe
;;;
;;; EXAMPLE:
;;; res = simulate_calib_source(300, 2000, 4000, 4500, 1.8E-05, 5.6E-5)  
;;;  
;;; MODIFICATION HISTORY:
;;;       Written by:     
;;;               Claudio Pagani

  
 if n_params() lt 6 then begin
    print,"Please enter correct parameters: res = simulate_calib_source(n_counts_in_line, detx, detymin, detymax, sctife, pctife)"
    print,"Please enter correct parameters: res = simulate_calib_source(300, 2000, 4000, 4500, 1.8E-05, 5.6E-5) EoL CTI estimates"
    print,"Please enter correct parameters: res = simulate_calib_source(300, 2000, 4000, 4500, 1.5E-06, 9.7E-7) BoL CTI estimates"
    return, "Try again please."
 endif  

  
  DEVICE,DECOMPOSED=0.
  tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1

  al_sigma = 30.0
  al_energy = 1486.0
  fe_energy = 5890.0

                                ; Generate uniformely distributed
                                ; dety coords in the segment of column
                                ; [detymin, detymax]
  
  dety = randomu(seed, n_counts_in_line)*(detymax-detymin)+detymin

                                ; Generate intrinsic energy line events, using a
                                ; standard deviation of 30eV at Al
                                ; K-alpha

   en_start = randomn(seed, n_counts_in_line)*al_sigma + al_energy
  
   sctial = sctife * (al_energy/fe_energy)^(-0.69623)
   pctial = pctife * (al_energy/fe_energy)^(-0.69623)
   
                                ; Calculate expected damaged energy
                                ; after transfers for each photon
   
   damaged_parallel = en_start * (1.-pctial)^dety
   damaged_serial = damaged_parallel * (1.-sctial)^detx

   losses_electrons_parallel = (en_start - damaged_parallel)/3.65
   losses_electrons_serial = (damaged_parallel - damaged_serial)/3.65
      
                                ; Randomise the process for each
                                ; photon, using a binomial process
                                ; with:
                                ; number of trials = readout
                                ; steps
                                ; Prob = probability of an electron capture in a pixel

   probability_capture_in_single_transfer_parallel = mean(losses_electrons_parallel/dety)
   probability_capture_in_single_transfer_serial = mean(losses_electrons_serial/detx)
   
   losses_electrons_stocastic_binomial = dblarr(n_counts_in_line)
   for count_photons = 0, n_counts_in_line-1 do losses_electrons_stocastic_binomial[count_photons] = randomu(seed, 1, binomial = [dety[count_photons], probability_capture_in_single_transfer_parallel]) + randomu(seed, 1, binomial = [detx, probability_capture_in_single_transfer_serial])
                                ; Check losses are not negative (can
                                ; happen if CTI values very low due to
                                ; the stocastic process)
   index_negative = where(losses_electrons_stocastic_binomial lt 0., n_negative)
   if n_negative gt 0 then losses_electrons_stocastic_binomial[index_negative] = 0.0
   
   damaged_stocastic_binomial = en_start - losses_electrons_stocastic_binomial * 3.65

   plothist, en_start, xhist_start, yhist_start, background=1, color=0, xr=[1000, 1600], bin = 5
   plothist, damaged_stocastic_binomial, xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, bin = 5, /overplot, color=3


                                ; Gaussian fit of input line
   res_gaussfit_start = gaussfit(xhist_start, yhist_start, coeff_fit_start, nterms=3, sigma = sigma_start)
   ploterror, xhist_start, yhist_start, sqrt(yhist_start), background=1, color=0, psym=1,  ERRCOLOR = 0
   oplot, xhist_start, res_gaussfit_start, color=3  

                                ; Gaussian fit of damaged line after transfers
   ploterror, xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, sqrt(yhist_damage_stocastic_bin), background=1, color=0, psym=1,  ERRCOLOR = 0
   res_gaussfit = gaussfit(xhist_damaged_stocastic_bin, yhist_damage_stocastic_bin, coeff_fit, nterms=3, sigma = sigma)
   oplot, xhist_damaged_stocastic_bin, res_gaussfit, color=3  
   cti_measured = 1.-(coeff_fit[1]/al_energy)^(1./(detx+(detymin+detymax)*0.5))

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
