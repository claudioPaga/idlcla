pro distort_cti_qdp_detx_dety_min_max_pixel_and_cti_correction, qdp_filename, expo_qdp, detx_min, detx_max, dety_min, dety_max, parallel_trap_density_pixel, serial_trap_density_pixel, cti_alpha, cti_beta, cti_corr_alpha, cti_corr_beta, serial_cti_alpha, serial_cti_beta, serial_cti_corr_alpha, serial_cti_corr_beta, cmax, threshold, key_release = key_release, key_plot = key_plot, detection_ll = detection_ll

;;;
;;; CP, 17 Dec 2018
;;; Summary - Distorts qdp file using a simple CTI model for a source
;;;           using a pixel-based implementation
;;;           Applies correction to distorted spectrum, using original
;;;           CTI parameters used to apply the distortion, that is,
;;;           CTI(Fe) = 10^-4 (SMILE EoL)
;;;           CTI beta = -0.25 (same as Swift)
;;;
;;; Input: qdp spectral file
;;;        expo_qdp - Exposure time used in XSPEC to simulate input
;;;        spectrum
;;;        dety_min, dety_max = source boundaries (Ex: 1250, 1400)
;;;        parallel_trap_density_pixel = number of traps in a pixel in
;;;        image section
;;;        serial_trap_density_pixel = number of traps in a pixel in
;;;        the readout register 
;;;        cti_alpha, cti_beta = parallel CTI distortion parameters
;;;        cti_corr_alpha, cti_corr_beta = CTI correction pars, not the
;;;        same as distortion, from SMILE EoL
;;;        preditions,
;;;        serial_cti_alpha, serial_cti_beta = serial CTI distortion parameters
;;;        serial_cti_corr_alpha, serial_cti_corr_beta = CTI correction pars, not the
;;;        same as distortion, from SMILE EoL
;;;        preditions,  
;;;        cmax, threshold = CTI correction iteration parameters, cmax
;;;        = max number of trials in iteration, threshold = min
;;;        differenct between derived value of corrected energy in the
;;;        iteration, in eV
;;;
;;; KEYWORD PARAMETERS:
;;;        key_release = recover energy within the 6-pixel binning,
;;;        default = 0
;;;        key_plot = generate plots
;;;        detection_ll = Min eV value to detect events (like lower  
;;;        level discriminator for XRT, it's the value to
;;;        "separate" real X-ray events from noise. Default = 0 (no threshold)
;;;  
;;;
;;; Output: qdp CTI distorted file/plots
;;;  
;;; HISTORY - Written by CP, 17 Dec 2018
;;;
;;;           Modified from distort_cti_qdp_dety_min_max_pixel
;;;           added procedure to correct for CTI losses after CTI
;;;           damage, to check how well energies can be recovered
;;;  
;;;           Modified from distort_cti_qdp_dety_min_max.pro
;;;           Distortion over a range of dety values, dety_min, dety_max
;;;           Distortion assumes observations are uniformely
;;;           distributed between dety_min and dety_max
;;;           Pixel-based treatment
;;;  
;;;           21/1/2019
;;;           Modified from distort_cti_qdp_dety_min_max_pixel_and_cti_correction.pro to implement serialCTI
;;;           distortion for a source imaged between detx_min and detx_max
;;;           and by added low threshold discriminator (keyword detection_ll) to
;;;           exclude low-energy X-rays that the CCD will not be able to
;;;           detect. it's the low level discriminator in XRT, to
;;;           separate noise from real X-rays  


  if n_params() lt 18 then begin
     print,"Please enter correct parameters: distort_cti_qdp_detx_dety_min_max_pixel_and_cti_correction, qdp_filename, expo_qdp, detx_min, detx_max, dety_min, dety_max, trap_density_pixel_parallel, trap_density_pixel_serial, cti_alpha, cti_beta, cti_corr_alpha, cti_corr_beta, serial_cti_alpha, serial_cti_beta, serial_cti_corr_alpha, serial_cti_corr_beta, cmax, threshold"
     print,"Please enter correct parameters: distort_cti_qdp_detx_dety_min_max_pixel_and_cti_correction, 'ccd2_sxi_puppis_strips_dety3893_expo25200.qdp', 25200, 10, 2250, 3893, 4029, 1., 1., 0.031, 0.44, 0.018843, -0.7, 0.012, 0.49, 0.00604860, -0.7, 10, 1, key_plot=1, detection_ll=100"
     return
  endif  

                                ; Trap release info
                                
  trap_species_fraction_parallel = [0.6, 0.4]
  trap_species_release_constant_parallel = [1., 1.]
  transfer_time_parallel = 1.
  trap_species_fraction_serial = [0.9, 0.1]
  trap_species_release_constant_serial = [1., 1.]
  transfer_time_serial = 1.
  if not(keyword_set(key_release)) then key_release = 0 ; Default = no recovery of released charge in pixel binning
  
  DEVICE,DECOMPOSED=0.
  tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1
  if not(keyword_set(detection_ll)) then detection_ll = 0.
  pixel_binning = 6
  readcol, qdp_filename, en, hw, counts_kev, err_kev, format='(d,d,d,d)', /silent

                                ;Derive best fit coeff of functional form of errors
                                ;This is needed for output qdp file.
                                ;(Funct form is a quadratic)
  

  res_lin = linfit(counts_kev, err_kev^2)
  counts_bin = round(counts_kev * hw * 2 * replicate(expo_qdp, n_elements(counts_kev)), /L64) ; replicate is needed in case the input parameter expo_qdp is a 1D array with 1 value (instead of a scalar). If replicate not used counts_bin will be a 1D array of 1 value, not what we want.
  counts_tot = long(total(counts_bin))
  list_energies = dblarr(counts_tot)
  list_energies_cti = dblarr(counts_tot)
  list_energies_cti_corrected = dblarr(counts_tot)
  losses_array_ev = dblarr(counts_tot)
  parallel_n_transfers_stored = dblarr(counts_tot)
  serial_n_transfers_stored = dblarr(counts_tot)
  parallel_cti_damaged_stored = dblarr(counts_tot)
  serial_cti_damaged_stored = dblarr(counts_tot)
  parallel_cti_corrections_ev_stored = dblarr(counts_tot)
  serial_cti_corrections_ev_stored = dblarr(counts_tot)
  parallel_recovered_ev_stored = dblarr(counts_tot)
  serial_recovered_ev_stored = dblarr(counts_tot)
  counter = 0l
  for i = 0, n_elements(en)-1 do begin
    for j = 0, counts_bin[i]-1 do begin
       list_energies[counter] = en[i]      
       ++counter
       ;print, counter
    endfor
 endfor

 ;;; Simple model of electron capture model.
 ;;; capture_prob = prob of capturing an electron from a charge packet
 ;;; of size list_energies
 ;;; losses_mean = expected electron losses, given the capture
 ;;; probability, the number of transfers and the density of traps in
 ;;; a pixel
 ;;; losses_sigma = standard deviation of expected caputes (model the
 ;;; probabilistic nature of trapping)

 ;;; Pixel-based treatment

  print, 'Number of X-rays in spectrum: ', counts_tot
  for count_events = 0, n_elements(list_energies)-1 do begin

    ; Parallel CTI
    eve_energy = list_energies[count_events]
    counter = 0
    ; Assign a DETY position to the X-ray event
    parallel_n_transfers = fix(randomu(Seed, 1) * (dety_max[0]-dety_min[0])) + dety_min[0]
    repeat begin
       ;Derive loss in pixel, update event energy
       capture_prob = 1. - exp(-cti_alpha * eve_energy^(1-cti_beta))
       ;if counter eq 0 then capture_prob_store = [capture_prob] else capture_prob_store = [capture_prob_store, capture_prob]
       ;Sample from binomial distribution
       res = randomn(seed, BINOMIAL = [parallel_trap_density_pixel , capture_prob])
       losses_kev = res * 3.65/1000.
       ;if losses_kev lt 0. then begin
       ;   print, capture_prob, losses_kev, format='(d15.12, d15.12)'
       ;   stop
       ;endif
       eve_energy -= losses_kev
       counter++
    endrep until (counter eq parallel_n_transfers or eve_energy le 0.)
    
                                ; Recovery of release charge in binning
    if key_release then begin
       pixels_released = parallel_n_transfers mod pixel_binning
       parallel_en_losses_kev = list_energies[count_events] - eve_energy 
       parallel_charge_released_within_binning =  total(trap_species_fraction_parallel  * parallel_en_losses_kev * (1. - exp(-pixels_released*transfer_time_parallel/ trap_species_release_constant_parallel)))
       parallel_recovered_ev_stored[count_events] = parallel_charge_released_within_binning*1000
       eve_energy += parallel_charge_released_within_binning    
    endif
    
    ;if eve_energy gt list_energies[count_events] then stop
    parallel_n_transfers_stored[count_events] = parallel_n_transfers
    parallel_cti_damaged_stored[count_events] = eve_energy

    ; Serial CTI
    serial_counter = 0
    ; Assign a DETX position to the X-ray event
    serial_n_transfers = fix(randomu(Seed, 1) * (detx_max[0]-detx_min[0])) + detx_min[0]
    serial_n_transfers_stored[count_events] = serial_n_transfers
    if eve_energy gt 0 then begin
       repeat begin
       ;Derive loss in pixel, update event energy
          serial_capture_prob = 1. - exp(-serial_cti_alpha * eve_energy^(1-serial_cti_beta))
          ;if serial_counter eq 0 then serial_capture_prob_store = [serial_capture_prob] else serial_capture_prob_store = [serial_capture_prob_store, serial_capture_prob]
                                ;Sample from binomial distribution
          res = randomn(seed, BINOMIAL = [serial_trap_density_pixel , serial_capture_prob])
          losses_kev = res * 3.65/1000.
          eve_energy -= losses_kev
          serial_counter++
       endrep until (serial_counter eq serial_n_transfers or eve_energy le 0.)

                                ; Recovery of release charge in binning
       if key_release then begin
          pixels_released = serial_n_transfers mod pixel_binning
          serial_en_losses_kev = parallel_cti_damaged_stored[count_events] - eve_energy 
          serial_charge_released_within_binning =  total(trap_species_fraction_serial  * serial_en_losses_kev * (1. - exp(-pixels_released*transfer_time_serial/ trap_species_release_constant_serial)))
          serial_recovered_ev_stored[count_events] = serial_charge_released_within_binning*1000
          eve_energy += serial_charge_released_within_binning    
       endif
       
                                ;plot, capture_prob_store
       list_energies_cti[count_events] = eve_energy
       ;serial_cti_losses_stored[count_events] = res * 3.65
    endif else begin
       list_energies_cti[count_events] = 0.
    endelse
    s_damaged = list_energies_cti[count_events]     
    serial_cti_damaged_stored[count_events] = list_energies_cti[count_events] 
    
    ; CTI Corrections

  ;  if eve_energy gt detection_ll/1000. then begin
    ; Serial correction       
  ;     eve_energy_ev = eve_energy * 1000.
                                ; Estimate event location (due to 6pixels binning)
 ;      x_location = floor((floor(serial_n_transfers/pixel_binning)+randomu(seed, 1))*pixel_binning+1.)
  ;     e_corr_serial_cti = cti_correction_function(eve_energy_ev, x_location, serial_cti_corr_alpha, serial_cti_corr_beta, cmax, threshold)
   ;    serial_cti_corrections_ev_stored[count_events] = e_corr_serial_cti - eve_energy_ev
  ;;                              ; Parallel correction
                                ; Estimate event location (due to 6pixels binning)
   ;    y_location = floor((floor(parallel_n_transfers/pixel_binning)+randomu(seed, 1))*pixel_binning+1.)
   ;    e_corr_parallel_cti = cti_correction_function(e_corr_serial_cti, y_location, cti_corr_alpha, cti_corr_beta, cmax, threshold)
   ;    list_energies_cti_corrected[count_events] = e_corr_parallel_cti/1000.
   ;    parallel_cti_corrections_ev_stored[count_events] = e_corr_parallel_cti - e_corr_serial_cti
       ;print, count_events, list_energies[count_events]*1000., p_damaged*1000., s_damaged*1000., e_corr_serial_cti, e_corr_parallel_cti, format='(i10, d10.4, d10.4, d10.4, d10.4, d10.4)'
       ;stop
   ; endif else begin
   ;    list_energies_cti_corrected[count_events] = 0.
   ;    serial_cti_corrections_ev_stored[count_events] = 0.
   ;    parallel_cti_corrections_ev_stored[count_events] = 0.
   ; endelse

 endfor

  
                                ; Derive serial and parallel CTI
                                ; correction values

 !P.Multi = [0, 1, 2, 0, 0]

  cti_parallel = 1. - (parallel_cti_damaged_stored/list_energies)^(1./parallel_n_transfers_stored)
  cti_serial = 1. - (serial_cti_damaged_stored/parallel_cti_damaged_stored)^(1./serial_n_transfers_stored)

                                ; Parallel CTI linear fit
                                ; Avoid CTI = 0 in fit, re-evaluate
                                ; CTI-correction parameters only if
                                ; there are enough events for a good
                                ; fit, otherwise use the input ones.
  index_values_ok_parallel = where(parallel_n_transfers_stored gt 100. and cti_parallel ne 0, n_fit_values_parallel)

  if n_fit_values_parallel gt 1000 then begin
  
     cti_parallel_logscale = alog10(cti_parallel[index_values_ok_parallel])
     list_energies_logscale_ev = alog10(list_energies[index_values_ok_parallel]*1000.)
     res_linear_parallel = linfit(list_energies_logscale_ev, cti_parallel_logscale, yfit = yfit_linfit_parallel)
     print, 'Parallel CTI fit: ', res_linear_parallel
     plot, list_energies_logscale_ev, cti_parallel_logscale, psym=1, color=0, background=1, ytit = 'Parallel CTI' 
     oplot, list_energies_logscale_ev, yfit_linfit_parallel, color=3

     print, 'Warning, new parameters for parallel CTI corrections were derived'
     print, 'Old values = ', cti_corr_alpha, cti_corr_beta
     cti_corr_alpha = 10^res_linear_parallel[0]
     cti_corr_beta = res_linear_parallel[1]
     print, 'New values = ', cti_corr_alpha, cti_corr_beta
         
  endif
     
     
                                ; Serial CTI linear fit
                                ; Avoid CTI = 0 in fit, re-evaluate
                                ; CTI-correction parameters only if
                                ; there are enough events for a good
                                ; fit, otherwise use the input ones.
  index_values_ok_serial = where(serial_n_transfers_stored gt 100. and cti_serial ne 0, n_fit_values_serial)
  if n_fit_values_serial gt 1000 then begin
   
     parallel_cti_damaged_stored_logscale_ev = alog10(parallel_cti_damaged_stored[index_values_ok_serial]*1000.)
     cti_serial_logscale = alog10(cti_serial[index_values_ok_serial])
     res_linear_serial = linfit(parallel_cti_damaged_stored_logscale_ev, cti_serial_logscale, yfit = yfit_linfit_serial)
     print, 'Serial CTI fit: ', res_linear_serial
     plot, parallel_cti_damaged_stored_logscale_ev, cti_serial_logscale, psym=1, color=0, background=1, ytit = 'Serial CTI'
     oplot, parallel_cti_damaged_stored_logscale_ev, yfit_linfit_serial, color=3

     print, 'Warning, new parameters for Serial CTI corrections were derived'
     print, 'Old values = ', serial_cti_corr_alpha, serial_cti_corr_beta
     serial_cti_corr_alpha = 10^res_linear_serial[0]
     serial_cti_corr_beta = res_linear_serial[1]
     print, 'New values = ', serial_cti_corr_alpha, serial_cti_corr_beta
         
    
     !P.Multi = [0, 1, 1, 0, 0]
  endif
  

    
; Remove events below the lower level discriminator, these events will be lost
 index_detection_ll_energies = where(list_energies_cti gt detection_ll/1000., n_ll) ;detection_ll is in eV, but list_energies_cti is in keV!
 print, 'N events, N above threshold after serial losses, total lost events: ', n_elements(list_energies_cti), n_ll,  n_elements(list_energies)-n_ll
 if n_ll gt 0 then begin
    list_energies_ll = list_energies[index_detection_ll_energies]
    list_energies_cti_ll = list_energies_cti[index_detection_ll_energies]
    parallel_n_transfers_stored_ll = parallel_n_transfers_stored[index_detection_ll_energies]
    serial_n_transfers_stored_ll = serial_n_transfers_stored[index_detection_ll_energies]
    parallel_cti_damaged_stored_ll = parallel_cti_damaged_stored[index_detection_ll_energies]
    serial_cti_damaged_stored_ll = serial_cti_damaged_stored[index_detection_ll_energies]
 endif else begin
    print, 'CTI very high, all events are lost below the low level discriminator!!!'
    return
 endelse
 
  cti_parallel_ll = 1. - (parallel_cti_damaged_stored_ll/list_energies_ll)^(1./parallel_n_transfers_stored_ll)
  cti_serial_ll = 1. - (serial_cti_damaged_stored_ll/parallel_cti_damaged_stored_ll)^(1./serial_n_transfers_stored_ll)
                                
  ; -----****---- START CTI Corrections
 ; A loop is needed as the corection function works on single energy events.
                                ; Estimate event location (due to 6pixels binning)
 list_energies_cti_corrected = dblarr(n_elements(list_energies_cti_ll))
 serial_cti_corrections_ev_stored = dblarr(n_elements(list_energies_cti_ll))
 parallel_cti_corrections_ev_stored = dblarr(n_elements(list_energies_cti_ll))
 
 for count_events = 0, n_elements(list_energies_cti_ll) - 1 do begin

   ; Serial correction       
       eve_energy_ev = list_energies_cti_ll[count_events] * 1000.
                                ; Estimate event location (due to 6pixels binning)
       x_location = floor((floor(serial_n_transfers_stored_ll[count_events]/pixel_binning)+randomu(seed, 1))*pixel_binning+1.)
       e_corr_serial_cti = cti_correction_function(eve_energy_ev, x_location, serial_cti_corr_alpha, serial_cti_corr_beta, cmax, threshold)
       serial_cti_corrections_ev_stored[count_events] = e_corr_serial_cti - eve_energy_ev
   ;                            ; Parallel correction
                                ; Estimate event location (due to 6pixels binning)
       y_location = floor((floor(parallel_n_transfers_stored[count_events]/pixel_binning)+randomu(seed, 1))*pixel_binning+1.)
       e_corr_parallel_cti = cti_correction_function(e_corr_serial_cti, y_location, cti_corr_alpha, cti_corr_beta, cmax, threshold)
       list_energies_cti_corrected[count_events] = e_corr_parallel_cti/1000.
       parallel_cti_corrections_ev_stored[count_events] = e_corr_parallel_cti - e_corr_serial_cti
       ;print, count_events, list_energies[count_events]*1000., p_damaged*1000., s_damaged*1000., e_corr_serial_cti, e_corr_parallel_cti, format='(i10, d10.4, d10.4, d10.4, d10.4, d10.4)'
       ;stop
     
 endfor

 ;---*******--- STOP CTI Corrections ;

 ;---*******--- PLOTS

 plothist, list_energies, xhist_input, yhist_input, bin=hw[0]*2, background=1, color=0
 plothist, list_energies_cti_ll, xhist, yhist, bin=hw[0]*2, /overplot, color=3
 plothist, list_energies_cti_corrected, xhist_corrected, yhist_corrected, bin=hw[0]*2, /overplot, color=4
 ;Determine ymax for postscipt plot
 ymaxplot = max([yhist_input, yhist, yhist_corrected])
 plothist, list_energies, xhist_input, yhist_input, bin=hw[0]*2, background=1, color=0, yr=[0, ymaxplot], boxplot=0
 plothist, list_energies_cti_ll, xhist, yhist, bin=hw[0]*2, /overplot, color=3, boxplot=0
 plothist, list_energies_cti_corrected, xhist_corrected, yhist_corrected, bin=hw[0]*2, /overplot, color=4, boxplot=0
 string_split = strsplit(qdp_filename, '.qdp', /regex, /extract)

                                ;Input CTI parameters are arbitrary,
                                ;if set too high all charge can be
                                ;lost. CTI qdp file and plot only
                                ;generated if there is some positive
                                ;chage (at least 10 bins)
 only_positive = where(xhist ge 0, n_positive)

 if n_positive gt 10 then begin

    xhist_gt0 = xhist[only_positive]
    yhist_gt0 = yhist[only_positive]

                                ;Estimate the errors in the
                                ;cti-distorted fluxes, based on the
                                ;best-fit parameters of the input
                                ;errors vs counts values
    err_cti_kev = sqrt(res_lin[0] + res_lin[1] * yhist_gt0)



    outfilenameqdp = 'pixel_cti_'+string_split[0]+'_detxmin'+strtrim(string(detx_min[0]),2)+'_detxmax'+strtrim(string(detx_max[0]),2)+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'_release'+strtrim(string(key_release),2)+'.qdp'
    openw, lu, outfilenameqdp, /get_lun
    print_array = [transpose(xhist_gt0), transpose(replicate(hw[0],n_elements(xhist_gt0))), transpose(yhist_gt0/(hw[0] * 2 * expo_qdp[0])), transpose(err_cti_kev)]
    printf, lu, print_array
    free_lun, lu

    ;Save cti-corrected spectrum
    xhist_corrected_gt0 = xhist_corrected[only_positive]
    yhist_corrected_gt0 = yhist_corrected[only_positive]
    err_cti_corrected_kev = sqrt(res_lin[0] + res_lin[1] * yhist_corrected_gt0)
    
    outfilenameqdp = 'pixel_cti_with_correction_'+string_split[0]+'_detxmin'+strtrim(string(detx_min[0]),2)+'_detxmax'+strtrim(string(detx_max[0]),2)+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'_release'+strtrim(string(key_release),2)+'.qdp'
    openw, lu, outfilenameqdp, /get_lun
    print_array = [transpose(xhist_corrected_gt0), transpose(replicate(hw[0],n_elements(xhist_corrected_gt0))), transpose(yhist_corrected_gt0/(hw[0] * 2 * expo_qdp[0])), transpose(err_cti_corrected_kev)]
    printf, lu, print_array
    free_lun, lu
    


    
    if keyword_set(key_plot) then begin
       !P.Multi = [0, 1, 1, 0, 0]
       entry_device = !d.name
       plotfilenameps = 'pixel_cti_with_correction_'+string_split[0]+'_detxmin'+strtrim(string(detx_min[0]),2)+'_detxmax'+strtrim(string(detx_max[0]),2)+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'_release'+strtrim(string(key_release),2)+'.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color
       titstring = 'Input [black] + distorted [red] + corrected [blue] '+ string_split[0]
       plothist, list_energies, xhist, yhist, bin=hw[0]*2, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "Cnts/s/keV", yr=[0, ymaxplot], boxplot=0, charsize=0.8,xr=[0.05,3]
       plothist, list_energies_cti_ll, bin=hw[0]*2, /overplot, color=3, boxplot=0
       plothist, list_energies_cti_corrected, bin=hw[0]*2, /overplot, color=4, boxplot=0
       expostr = 'Expo = '+ strtrim(string(floor(expo_qdp[0])), 2)+ ' s'
       xyouts, 0.7, 0.8, expostr, /norm, charsize=0.8
       driftstr = 'Drift = '+ strtrim(string(floor(dety_max[0]-dety_min[0])), 2)+ ' pixels'
       xyouts, 0.7, 0.75, driftstr, /norm, charsize=0.8
       xtrmaxstr = 'DETX = '+ strtrim(string(floor(detx_min[0])), 2)+ '-' +  strtrim(string(floor(detx_max[0])), 2)
       xyouts, 0.7, 0.7, xtrmaxstr, /norm, charsize=0.8
       ytrmaxstr = 'DETY = '+ strtrim(string(floor(dety_min[0])), 2)+ '-' +  strtrim(string(floor(dety_max[0])), 2)
       xyouts, 0.7, 0.65, ytrmaxstr, /norm, charsize=0.8
       device,/close  
       set_plot,entry_device

       
       plotfilenameps = 'pixel_cti_with_correction_'+string_split[0]+'_detxmin'+strtrim(string(detx_min[0]),2)+'_detxmax'+strtrim(string(detx_max[0]),2)+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'_release'+strtrim(string(key_release),2)+'_log.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color
       titstring = 'Input [black] + distorted [red] + corrected [blue] ' + string_split[0]
       plothist, list_energies, bin=hw[0]*2, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "Cnts/s/keV",xr=[0.05,3],/ylog, yr=[1, ymaxplot], xst=1, yst=1, boxplot=0, charsize=0.8
       plothist, list_energies_cti_ll, bin=hw[0]*2, /overplot, color=3, boxplot=0
       plothist, list_energies_cti_corrected, bin=hw[0]*2, /overplot, color=4, boxplot=0
       xyouts, 0.7, 0.8, expostr, /norm, charsize=0.8
       xyouts, 0.7, 0.75, driftstr, /norm, charsize=0.8
       xyouts, 0.7, 0.7, xtrmaxstr, /norm , charsize=0.8
       xyouts, 0.7, 0.65, ytrmaxstr, /norm, charsize=0.8
       device,/close  
       set_plot,entry_device

       ;Additional plots of distribution of losses vs DETY and ENERGY
       ;!P.Multi = [0, 1, 2, 0, 0]
       ;plotfilenameps = 'pixel_losses_vs_dety_energy_'+string_split[0]+'.ps'
       ;set_plot,'ps'
       ;device, filename = plotfilenameps, /color
       ;titstring = 'CTI losses '+ string_split[0]
       ;plot, list_energies, losses_array_ev, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "CTI Losses (eV)",xr=[0.05,3], psym=1
       ;plot, serial_n_transfers, losses_array_ev, background=1, color=0, tit = titstring, xtit="DETY", ytit = "CTI Losses (eV)", psym=1       
       ;device,/close  
       ;set_plot,entry_device
       ;!P.Multi = [0, 1, 1, 0, 0]
   
    endif

 endif else print, 'WARNING - CTI TOO HIGH, ALL CHARGE IS LOST!!!!!'

 ;en6kev = 6.
 ;losses6keV= (1. - exp(-cti_alpha * en6kev^(1-cti_beta))) * n_transfers * trap_density_pixel * 3.65
 ;en1kev = 1.
 ;losses1keV= (1. - exp(-cti_alpha * en1kev^(1-cti_beta))) * n_transfers * trap_density_pixel * 3.65
 ;en500ev = 0.5
 ;losses500ev = (1. - exp(-cti_alpha * en500ev^(1-cti_beta))) * n_transfers * trap_density_pixel * 3.65
 
 ;losses_corrections = losses_array_ev - cti_corrections_ev_stored
 ;plot, list_energies, losses_corrections, psym=1, background=1, color=0, ytit='Losses-corrections'
 ;print, 'Mean and median of losses-corrections'
 ;print, mean(losses_corrections, /nan), median(losses_corrections)
  
; print, 'Losses at 6 keV = ' ,losses6keV
; print, 'Losses at 1 keV = ' ,losses1keV
; print, 'Losses at .5 keV = ',losses500ev
end

 









