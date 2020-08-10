pro distort_cti_qdp_and_cti_correction, qdp_filename, expo_qdp, n_transfers, trap_density_pixel, cti_alpha, cti_beta, cti_corr_alpha, cti_corr_beta, cmax, threshold, key_plot = key_plot, key_serial = key_serial, serial_n_transfers = serial_n_transfers, serial_cti_alpha = serial_cti_alpha, serial_cti_beta = serial_cti_beta, serial_cti_corr_alpha = serial_cti_corr_alpha, serial_cti_corr_beta = serial_cti_corr_beta, detection_ll = detection_ll 

;;;
;;; CP, 18 Dec 2018
;;; Summary - Distorts qdp file using a simple CTI model
;;;           Applies correction to distorted spectrum, using original
;;;           CTI parameters used to apply the distortion, from OU
;;;           estimate;  
;;;           pCTI(E) = 0.018843e^-0.7
;;;           This translates into a CTI at Fe of  5.610^-5 (OU SMILE EoL)
;;;           sCTI(E) = 0.321 pCTI(E)
;;;  
;;; Input: qdp spectral file
;;;        expo_qdp - Exposure time used in XSPEC to simulate input
;;;        spectrum
;;;        n_transfers = source average DETY position
;;;        trap_density_pixel = number of traps in a pixel
;;;        cti_alpha, cti_beta = CTI distortion parameters
;;;        cti_corr_alpha, cti_corr_beta = CTI correction pars, not the
;;;        same as distortion, these should be 0.018843, 0.7 as in SMILE EoL
;;;        preditions
;;;        cmax, threshold = CTI correction iteration parameters, cmax
;;;        = max number of trials in iteration, threshold = min
;;;        differenct between derived value of corrected energy in the
;;;        iteration, in eV
;;; KEYWORD PARAMETERS:
;;;        key_plot = generate plots
;;;        key_serial = flag to implements/not serial CTI distortion. If not
;;;        provided default is key_serial = 0
;;;        serial_n_transfers = DETX position (should be between 1 and 2255
;;;        (because 2 output notdes)), default = 1000
;;;        serial_cti_alpha = 0.012
;;;        serial_cti_beta =  0.49
;;;        serial_cti_corr_alpha = 0.321*0.018843
;;;        serial_cti_corr_beta = 0.7 
;;;        detection_ll = Min eV value to detect events (like lower
;;;        level discriminator for XRT, it's the value to
;;;        "separate" real X-ray events from noise. Default = 0 (no threshold)
;;;  
;;; Output: qdp CTI distorted file
;;;
;;; HISTORY - Written by CP, 18 Dec 2018
;;;
;;;           Modified from distort_cti_qdp
;;;           added procedure to correct for CTI losses after CTI
;;;           damage, to check how well energies can be recovered
;;;  
;;; CP, 17/1/2019
;;; V2, Modified to implement serialCTI distortion
;;;
;;; CP, 18/1/2019
;;; V3, added low threshold discriminator (keyword detection_ll) to
;;; exclude low-energy X-rays that the CCD will not be able to
;;; detect. it's the low level discriminator in XRT, to
;;; separate noise from real X-rays  
  
 if n_params() lt 10 then begin
    print,"Please enter correct parameters: distort_cti_qdp_and_cti_correction, qdp_filename, expo_qdp, n_transfers, trap_density_pixel, cti_alpha, cti_beta, cti_corr_alpha, cti_corr_beta, cmax, threshold"
    print,"Please enter correct parameters: distort_cti_qdp_and_cti_correction, 'smile_sxi_casa_10ks_323.qdp', 10000., 1550, 1., 0.031, 0.44, 0.018843, -0.7, 10, 1, key_plot = 1, key_serial = 1, serial_n_transfers = 1000, serial_cti_alpha = 0.012, serial_cti_beta = 0.49, serial_cti_corr_alpha = 0.00604860, serial_cti_corr_beta = -0.7, detection_ll=100 "
    return
 endif  

  DEVICE,DECOMPOSED=0.
  tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1
  pixel_binning = 6

   ; set default serial CTI parameters
  if not(keyword_set(key_serial)) then key_serial = 0
  if not(keyword_set(serial_n_transfers)) then serial_n_transfers = 1000
  if not(keyword_set(serial_cti_alpha)) then serial_cti_alpha = 0.012
  if not(keyword_set(serial_cti_beta)) then serial_cti_beta = 0.49
  if not(keyword_set(serial_cti_corr_alpha)) then serial_cti_corr_alpha = 0.00604860
  if not(keyword_set(serial_cti_corr_beta)) then serial_cti_corr_beta = 0.7
  if not(keyword_set(detection_ll)) then detection_ll = 0.

  
  readcol, qdp_filename, en, hw, counts_kev, err_kev, format='(d,d,d,d)'

                                ;Derive best fit coeff of functional form of errors
                                ;This is needed for output qdp file.
                                ;(Funct form is a quadratic)
  

  res_lin = linfit(counts_kev, err_kev^2)
  counts_bin = round(counts_kev * hw * 2 * expo_qdp)
  counts_tot = long(total(counts_bin))
  list_energies = dblarr(counts_tot)
  list_energies_cti = dblarr(counts_tot)
  losses_array_ev = dblarr(counts_tot)
  list_energies_cti_corrected = dblarr(counts_tot)
  cti_losses_ev_storage = dblarr(counts_tot)
  

  
  counter = 0l
 for i = 0, n_elements(en)-1 do begin
    for j = 0, counts_bin[i]-1 do begin
       list_energies[counter] = en[i]      
       ++counter
;       print, counter
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
 
 capture_prob = 1. - exp(-cti_alpha * list_energies^(1-cti_beta))
 losses_mean = n_transfers * capture_prob * trap_density_pixel 
 losses_sigma = sqrt(n_transfers * trap_density_pixel * capture_prob * (1. - capture_prob))
 random_array = randomn(seed, n_elements(list_energies))
 losses_array_ev = (random_array * losses_sigma + losses_mean) * 3.65
                                  ;Avoid "negative" losses (can be
                                ;caused by low losses_mean value and a
                                ;large losses_sigma value
 index_l0 = where(losses_array_ev lt 0., nl0)
 if nl0 gt 0 then losses_array_ev[index_l0] = 0
 
 list_energies_cti  = (list_energies*1000. - losses_array_ev)/1000. ; Damaged flux in eV
 
; Remove negative energies, these events will be lost
 index_positive_energies = where(list_energies_cti gt 0., n_positive)
 if n_positive gt 0 then list_energies_cti = list_energies_cti[index_positive_energies] else begin
    print, 'CTI very high, all events are lost!!!'
    return
 endelse


 ;;; Apply serial CTI distortion (to distorted energies list_energies_cti

 if key_serial then begin

    serial_capture_prob = 1. - exp(-serial_cti_alpha * list_energies_cti^(1-serial_cti_beta))
    serial_losses_mean = serial_n_transfers * serial_capture_prob * trap_density_pixel 
    serial_losses_sigma = sqrt(serial_n_transfers * trap_density_pixel * serial_capture_prob * (1. - serial_capture_prob))
    serial_random_array = randomn(seed, n_elements(list_energies_cti))
    serial_losses_array_ev = (serial_random_array * serial_losses_sigma + serial_losses_mean) * 3.65
                                  ;Avoid "negative" losses (can be
                                ;caused by low losses_mean value and a
                                ;large losses_sigma value
    index_l0 = where(serial_losses_array_ev lt 0., nl0)
    if nl0 gt 0 then serial_losses_array_ev[index_l0] = 0
 
 ;;; Apply serial losses, avoid negative energies
    list_energies_cti  = (list_energies_cti*1000. - serial_losses_array_ev)/1000.
    index_negative_energies = where(list_energies_cti lt 0., n_negative)
    if n_negative gt 0 then list_energies_cti[index_negative_energies] = 0.
 endif

; Remove events below the lower level discriminatornegative energies, these events will be lost
 index_detection_ll_energies = where(list_energies_cti gt detection_ll/1000., n_ll) ;detection_ll is in eV, but list_energies_cti is in keV!
 if n_ll gt 0 then list_energies_cti = list_energies_cti[index_detection_ll_energies]  else begin
    print, 'CTI very high, all events are lost below the low level discriminator!!!'
    return
 endelse


 
 
 ; CTI Corrections
 ; A loop is needed as the corection function works on single energy events.
 ; Estimate event location (due to 6pixels binning)
 cti_corrections_ev_storage = dblarr(n_elements(list_energies_cti))
 print, 'N. events in spectrum: ', n_elements(list_energies_cti)
 for count_events = 0, n_elements(list_energies_cti) - 1 do begin
 ; Serial correction first
    eve_energy_ev = list_energies_cti[count_events] * 1000.
    e_corr_serial_cti = cti_correction_function(eve_energy_ev, serial_n_transfers, serial_cti_corr_alpha, serial_cti_corr_beta, cmax, threshold)
 ; Parallel correction
    location = floor((floor(n_transfers/pixel_binning)+randomu(seed, 1))*pixel_binning+1.)
    e_corr_cti = cti_correction_function(e_corr_serial_cti, location, cti_corr_alpha, cti_corr_beta, cmax, threshold)
    list_energies_cti_corrected[count_events] = e_corr_cti/1000.
    cti_corrections_ev_storage[count_events] = e_corr_cti - eve_energy_ev
                                ;if
                                ;losses_array_ev[count_events]-cti_corrections_ev_storage[count_events]
                                ;gt 30. then stop
    ;print, count_events, eve_energy_ev, e_corr_serial_cti, e_corr_cti
    
 endfor
 plothist, list_energies, xhist_input, yhist_input, bin=hw[0]*2, background=1, color=0, boxplot=0
 plothist, list_energies_cti, xhist, yhist, bin=hw[0]*2, /overplot, color=3, boxplot=0
 plothist, list_energies_cti_corrected, xhist_corrected, yhist_corrected, bin=hw[0]*2, /overplot, color=4, boxplot=0

                                ;Determine ymax for postscipt plot
 ymaxplot = max([yhist_input, yhist, yhist_corrected])
 plothist, list_energies, xhist_input, yhist_input, bin=hw[0]*2, background=1, color=0, yr=[0, ymaxplot], boxplot=0
 plothist, list_energies_cti, xhist, yhist, bin=hw[0]*2, /overplot, color=3, boxplot=0
 plothist, list_energies_cti_corrected, xhist_corrected, yhist_corrected, bin=hw[0]*2, /overplot, color=4, boxplot=0

 
 
 string_split = strsplit(qdp_filename, '.qdp', /regex, /extract)
 outfile_name = string_split[0]+'_cti_dety'+strtrim(string(n_transfers),2)+'.qdp'

                                ;Input CTI parameters are arbitrary,
                                ;if set too high all charge can be
                                ;lost. CTI qdp file and plot only
                                ;generated if there is some positive
                                ;chage (at least 10 bins)
 only_positive = where(xhist ge 0., n_positive)

 if n_positive gt 10 then begin

    xhist_gt0 = xhist[only_positive]
    yhist_gt0 = yhist[only_positive]

                                ;Estimate the errors in the
                                ;cti-distorted fluxes, based on the
                                ;best-fit parameters of the input
                                ;errors vs counts values
    err_cti_kev = sqrt(res_lin[0] + res_lin[1] * yhist_gt0)
 
    if key_serial then outfilenameqdp = 'cti_'+string_split[0]+'_detx'+strtrim(string(serial_n_transfers),2)+'_dety'+strtrim(string(n_transfers),2)+'.qdp' else outfilenameqdp = 'cti_'+string_split[0]+'_dety'+strtrim(string(n_transfers),2)+'.qdp'
    openw, lu, outfilenameqdp, /get_lun
    printf, lu, [transpose(xhist_gt0), transpose(replicate(hw[0],n_elements(xhist_gt0))), transpose(yhist_gt0/(hw * 2 * expo_qdp)), transpose(err_cti_kev)]
    free_lun, lu

;Save cti-corrected spectrum
    xhist_corrected_gt0 = xhist_corrected[only_positive]
    yhist_corrected_gt0 = yhist_corrected[only_positive]
    err_cti_corrected_kev = sqrt(res_lin[0] + res_lin[1] * yhist_corrected_gt0)
    if key_serial then outfilenameqdp = 'cti_with_correction_'+string_split[0]+'_detx'+strtrim(string(serial_n_transfers),2)+'_dety'+strtrim(string(n_transfers),2)+'.qdp' else outfilenameqdp = 'cti_with_correction_'+string_split[0]+'_dety'+strtrim(string(n_transfers),2)+'.qdp'
    openw, lu, outfilenameqdp, /get_lun
    print_array = [transpose(xhist_corrected_gt0), transpose(replicate(hw[0],n_elements(xhist_corrected_gt0))), transpose(yhist_corrected_gt0/(hw[0] * 2 * expo_qdp[0])), transpose(err_cti_corrected_kev)]
    printf, lu, print_array
    free_lun, lu
    
    if keyword_set(key_plot) then begin
       entry_device = !d.name
       if key_serial then plotfilenameps = 'cti_with_correction_'+string_split[0]+'_detx'+strtrim(string(serial_n_transfers),2)+'_dety'+strtrim(string(n_transfers),2)+'.ps'  else plotfilenameps = 'cti_with_correction_'+string_split[0]+'_dety'+strtrim(string(n_transfers),2)+'.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color,/tt_font, set_font='Times', font_size=10
       titstring = 'Input [black] + distorted [red] + corrected [blue] '+ string_split[0]
       plothist, list_energies, xhist, yhist, bin=hw[0]*2, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "Cnts/s/keV", yr=[0, ymaxplot], boxplot=0, charsize=0.8
       plothist, list_energies_cti, bin=hw[0]*2, /overplot, color=3, boxplot=0
       plothist, list_energies_cti_corrected, bin=hw[0]*2, /overplot, color=4, boxplot=0
       expostr = 'Expo = '+ strtrim(string(floor(expo_qdp[0])), 2)+ ' s'
       xyouts, 0.75, 0.8, expostr, /norm, charsize=1.;1
       trmaxstr = 'DETY = '+ strtrim(string(floor(n_transfers)), 2)
       xyouts, 0.75, 0.75, trmaxstr, /norm, charsize=1.;1
              
       device,/close  
       set_plot,entry_device
       
       if key_serial then plotfilenameps = 'cti_with_correction_'+string_split[0]+'_detx'+strtrim(string(serial_n_transfers),2)+'_dety'+strtrim(string(n_transfers),2)+'_log.ps'  else plotfilenameps = 'cti_with_correction_'+string_split[0]+'_dety'+strtrim(string(n_transfers),2)+'_log.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color
       titstring = 'Input [black] + distorted [red] + corrected [blue] ' + string_split[0]
       plothist, list_energies, bin=hw[0]*2, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "Cnts/s/keV",xr=[0.2,3],/ylog, yr=[1, ymaxplot], xst=1, yst=1, boxplot=0, charsize=0.8
       plothist, list_energies_cti, bin=hw[0]*2, /overplot, color=3, boxplot=0
       plothist, list_energies_cti_corrected, bin=hw[0]*2, /overplot, color=4, boxplot=0
       xyouts, 0.75, 0.8, expostr, /norm, charsize=1.;1
       xyouts, 0.75, 0.75, trmaxstr, /norm, charsize=1.;1
       device,/close  
       set_plot,entry_device

                                ;Side by side linear version

       !P.Multi = [0, 2, 2, 0, 0]
       entry_device = !d.name
       if key_serial then plotfilenameps = 'cti_with_correction_'+string_split[0]+'_detx'+strtrim(string(serial_n_transfers),2)+'_dety'+strtrim(string(n_transfers),2)+'_split.ps'  else plotfilenameps = 'cti_with_correction_'+string_split[0]+'_dety'+strtrim(string(n_transfers),2)+'_split.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color,/tt_font, set_font='Times', font_size=10
       plothist, list_energies, xhist, yhist, bin=hw[0]*2, background=1, color=0, xtit="Energy [keV]", ytit = "Input Spectrum (Cnts/s/keV)", boxplot=0
       plothist, list_energies_cti, bin=hw[0]*2, background=1, color=0, xtit="Energy [keV]", ytit = "CTI Spectrum (Cnts/s/keV)", boxplot=0
       plothist, list_energies_cti, bin=hw[0]*2, /overplot, color=3, boxplot=0
       plothist, list_energies_cti_corrected, bin=hw[0]*2, /overplot, color=4, boxplot=0
       
       device,/close  
       set_plot,entry_device

       !P.Multi = [0, 1, 1, 0, 0]
       
    endif

 endif else print, 'WARNING - CTI TOO HIGH, ALL CHARGE IS LOST!!!!!'

 en6kev = 6.
 losses6keV= (1. - exp(-cti_alpha * en6kev^(1-cti_beta))) * n_transfers * trap_density_pixel * 3.65
 en1kev = 1.
 losses1keV= (1. - exp(-cti_alpha * en1kev^(1-cti_beta))) * n_transfers * trap_density_pixel * 3.65
 en500ev = 0.5
 losses500ev = (1. - exp(-cti_alpha * en500ev^(1-cti_beta))) * n_transfers * trap_density_pixel * 3.65

 losses_corrections = losses_array_ev - cti_corrections_ev_storage
 plot, list_energies, losses_corrections, psym=1, background=1, color=0, ytit='Losses-corrections'
 print, 'Mean and median of losses-corrections'
 print, mean(losses_corrections, /nan), median(losses_corrections)
 stop
; print, 'Losses at 6 keV = ' ,losses6keV
; print, 'Losses at 1 keV = ' ,losses1keV
; print, 'Losses at .5 keV = ',losses500ev
end

 
