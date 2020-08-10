pro distort_cti_qdp_dety_min_max, qdp_filename, expo_qdp, dety_min, dety_max, trap_density_pixel, cti_alpha, cti_beta, key_plot = key_plot, key_serial = key_serial, serial_n_transfers = serial_n_transfers, serial_cti_alpha = serial_cti_alpha, serial_cti_beta = serial_cti_beta
;;;
;;;
;;; CP, 3 Dec 2018
;;; Summary - Distorts qdp file using a simple CTI model
;;;  
;;; Input: qdp spectral file
;;; Output: qdp CTI distorted file
;;;
;;; KEYWORD PARAMETERS:
;;; key_plot = generate plots
;;; key_serial = flag to implements/not serial CTI distortion. If not
;;; provided default is key_serial = 0
;;; serial_n_transfers = DETX position (should be between 1 and 2255
;;; (because 2 output notdes)), default = 1000
;;; serial_cti_alpha = 0.012 by default
;;; serial_cti_beta =  0.049 by default
;;;
;;;  
;;; HISTORY - Written by CP, 3 Dec 2018
;;;           Modified from distort_cti_qdp_pixel.pro
;;;           Distortion over a range of dety values, dety_min, dety_max
;;;           Distortion assumes observations are uniformely
;;;           distributed between dety_min and dety_max
;;; 
;;; V2, Modified to implement serialCTI distortion

  


  
 if n_params() lt 6 then begin
    print,"Please enter correct parameters: distort_cti_qdp_dety_min_max, qdp_filename, expo_qdp, dety_min, dety_max, trap_density_pixel, cti_alpha, cti_beta"
    print,"Please enter correct parameters: distort_cti_qdp_dety_min_max, 'smile_sxi_casa_10ks_323.qdp', 10000, 1000, 1550, 1., 0.1, 0.5"
    return
 endif  

  DEVICE,DECOMPOSED=0.
  tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1

  ; set default serial CTI parameters
  if not(keyword_set(key_serial)) then key_serial = 0
  if not(keyword_set(serial_n_transfers)) then serial_n_transfers = 1000
  if not(keyword_set(serial_cti_alpha)) then serial_cti_alpha = 0.012
  if not(keyword_set(serial_cti_beta)) then serial_cti_beta = 0.49

  
  readcol, qdp_filename, en, hw, counts_kev, err_kev, format='(d,d,d,d)'

                                ;Derive best fit coeff of functional form of errors
                                ;This is needed for output qdp file.
                                ;(Funct form is a quadratic)
  

  res_lin = linfit(counts_kev, err_kev^2)
  counts_bin = round(counts_kev * hw * 2 * replicate(expo_qdp, n_elements(counts_kev)), /L64) ; replicate is needed in case the input parameter expo_qdp is a 1D array with 1 value (instead of a scalar). If replicate not used counts_bin will be a 1D array of 1 value, not what we want.
  counts_tot = long(total(counts_bin))
  list_energies = dblarr(counts_tot)
  list_energies_cti = dblarr(counts_tot)
  losses_array_ev = dblarr(counts_tot)
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
 
 capture_prob = 1. - exp(-cti_alpha * list_energies^(1-cti_beta))

                                ;Distribution of number of transfers
                                ;(this corresponds to the DETY
                                ;position of the X-ray event. The
                                ;assumption is of a uniform
                                ;distribution along the column.
 n_transfers = fix(randomu(Seed, counts_tot) * (dety_max[0]-dety_min[0])) + dety_min[0]
 losses_mean = n_transfers * capture_prob * trap_density_pixel 
 losses_sigma = sqrt(n_transfers * trap_density_pixel * capture_prob * (1. - capture_prob))
 random_array = randomn(seed, n_elements(list_energies)-1)
 losses_array_ev = (random_array * losses_sigma + losses_mean) * 3.65 

                                ;Avoid "negative" losses (can be
                                ;caused by low losses_mean value and a
                                ;large losses_sigma value
 index_l0 = where(losses_array_ev lt 0., nl0)
 if nl0 gt 0 then losses_array_ev[index_l0] = 0
 

 list_energies_cti  = (list_energies*1000. - losses_array_ev)/1000. ; Damaged flux in keV

; Remove negative energies, these events will be lost
 index_positive_energies = where(list_energies_cti lt 0., n_positive)
 if n_positive gt 0 then list_energies_cti[index_positive_energies] = 0. else begin
    print, 'CTI very high, all events are lost!!!'
    return
 endelse


 ;;; Apply serial CTI distortion (to distorted energies list_energies_cti

 if key_serial eq 1 then begin
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

; orig = histogram(list_energies, bin = hw[0] * 2)
; cti = histogram(list_energies_cti, bin = hw[0] * 2, locations = xbin_cti_start)
 
; plot, en , counts_kev, BACKGROUND=1,color=0
; oplot, xbin_cti_start, cti, linestyle = 3, color = 4, thick=2


 plothist, list_energies, xhist_input, yhist_input, bin=hw[0]*2, background=1, color=0, box=0
 plothist, list_energies_cti, xhist, yhist, bin=hw[0]*2, /overplot, color=4, box=0
 
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

    if key_serial then  outfilenameqdp = 'cti_'+string_split[0]+'_detx'+strtrim(string(serial_n_transfers),2)+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'.qdp' else outfilenameqdp = 'cti_'+string_split[0]+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'.qdp'
    openw, lu, outfilenameqdp, /get_lun
    print_array = [transpose(xhist_gt0), transpose(replicate(hw[0],n_elements(xhist_gt0))), transpose(yhist_gt0/(hw[0] * 2 * expo_qdp[0])), transpose(err_cti_kev)]
    printf, lu, print_array
    free_lun, lu

    if keyword_set(key_plot) then begin
       entry_device = !d.name

       if key_serial then plotfilenameps = 'cti_'+string_split[0]+'_detx'+strtrim(string(serial_n_transfers),2)+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'.ps' else plotfilenameps = 'cti_'+string_split[0]+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'.ps'
       
       set_plot,'ps'
       device, filename = plotfilenameps, /color,/tt_font, set_font='Times', font_size=10
       titstring = 'Input [black] + distorted [red] '+ string_split[0]
       plothist, list_energies, xhist, yhist, bin=hw[0]*2, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "Cnts/s/keV", box=0, charsize=0.8
       plothist, list_energies_cti, bin=hw[0]*2, /overplot, color=3, box=0
       expostr = 'Expo = '+ strtrim(string(floor(expo_qdp[0])), 2)+ ' s'
       xyouts, 0.75, 0.8, expostr, /norm, charsize=1;.2
       driftstr = 'Drift = '+ strtrim(string(floor(dety_max[0]-dety_min[0])), 2)+ ' pixels'
       xyouts, 0.75, 0.75, driftstr, /norm, charsize=1;.2
       trmaxstr = 'DETY max = '+ strtrim(string(floor(dety_max[0])), 2)
       xyouts, 0.75, 0.7, trmaxstr, /norm, charsize=1;.2
       
       device,/close  
       set_plot,entry_device

       if key_serial then plotfilenameps = 'cti_'+string_split[0]+'_detx'+strtrim(string(serial_n_transfers),2)+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'_log.ps'  else plotfilenameps = 'cti_'+string_split[0]+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'_log.ps'
             
       set_plot,'ps'
       device, filename = plotfilenameps, /color
       titstring = 'Input [black] + distorted [red] ' + string_split[0]
       plothist, list_energies, bin=hw[0]*2, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "Cnts/s/keV",xr=[0.2,3],/ylog, yr=[1, max(yhist)*1.1], xst=1, yst=1, box=0, charsize=0.8
       plothist, list_energies_cti, bin=hw[0]*2, /overplot, color=3, box=0
       xyouts, 0.75, 0.8, expostr, /norm, charsize=1;.2
       xyouts, 0.75, 0.75, driftstr, /norm, charsize=1;.2
       xyouts, 0.75, 0.7, trmaxstr, /norm, charsize=1;.2
       device,/close  
       set_plot,entry_device

       ;Additional plots of distribution of losses vs DETY and ENERGY

       !P.Multi = [0, 1, 2, 0, 0]
        
       plotfilenameps = 'losses_vs_dety_energy_'+string_split[0]+'.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color
       titstring = 'CTI losses '+ string_split[0]
       plot, list_energies, losses_array_ev, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "CTI Losses (eV)",xr=[0.2,3], psym=1, charsize=0.8
       plot, n_transfers, losses_array_ev, background=1, color=0, tit = titstring, xtit="DETY", ytit = "CTI Losses (eV)", psym=1
       
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

; print, 'Losses at 6 keV = ' ,losses6keV
; print, 'Losses at 1 keV = ' ,losses1keV
; print, 'Losses at .5 keV = ',losses500ev
end

 
