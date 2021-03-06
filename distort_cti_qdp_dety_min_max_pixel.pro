pro distort_cti_qdp_dety_min_max_pixel, qdp_filename, expo_qdp, dety_min, dety_max, trap_density_pixel, cti_alpha, cti_beta, key_plot = key_plot
;;;
;;; CP, 6 Dec 2018
;;; Summary - Distorts qdp file using a simple CTI model for a source
;;;           observed within dety_min and dety_max, using a
;;;           pixel-based implementation
;;;  
;;; Input: qdp spectral file
;;; Output: qdp CTI distorted file
;;;  
;;; HISTORY - Written by CP, 6 Dec 2018
;;;           Modified from distort_cti_qdp_dety_min_max.pro
;;;           Distortion over a range of dety values, dety_min, dety_max
;;;           Distortion assumes observations are uniformely
;;;           distributed between dety_min and dety_max
;;;           Pixel-based treatment  



  
 if n_params() lt 6 then begin
    print,"Please enter correct parameters: distort_cti_qdp_dety_min_max_piexl, qdp_filename, expo_qdp, dety_min, dety_max, trap_density_pixel, cti_alpha, cti_beta"
    print,"Please enter correct parameters: distort_cti_qdp_dety_min_max_pixel, 'smile_sxi_casa_10ks_323.qdp', 10000, 1000, 1550, 1., 0.1, 0.5"
    return
 endif  

  DEVICE,DECOMPOSED=0.
  tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1


  readcol, qdp_filename, en, hw, counts_kev, err_kev, format='(d,d,d,d)', /silent

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

 ;;; Pixel-based treatment

  print, 'Number of X-rays in spectrum: ', counts_tot
  for count_events = 0, n_elements(list_energies)-1 do begin
    eve_energy = list_energies[count_events]
    counter = 0
    ; Assign a DETY position to the X-ray event
    n_transfers = fix(randomu(Seed, 1) * (dety_max[0]-dety_min[0])) + dety_min[0]
    repeat begin
       ;Derive loss in pixel, update event energy
       capture_prob = 1. - exp(-cti_alpha * eve_energy^(1-cti_beta))
       if counter eq 0 then capture_prob_store = [capture_prob] else capture_prob_store = [capture_prob_store, capture_prob]
       ;Sample from binomial distribution
       res = randomn(seed, BINOMIAL = [trap_density_pixel , capture_prob])
       losses_kev = res * 3.65/1000.
       eve_energy -= losses_kev
       counter++
    endrep until (counter eq n_transfers or eve_energy le 0.)
    ;plot, capture_prob_store
    list_energies_cti[count_events] = eve_energy
 endfor
       
 lost = where(list_energies_cti le 0., xrays_lost_count)
 print, 'Number of events losing all energy: ', xrays_lost_count
 
                                ;Avoid "negative" losses (can be
                                ;caused by low losses_mean value and a
                                ;large losses_sigma value
 index_l0 = where(list_energies_cti lt 0., nl0)
 if nl0 gt 0 then list_energies_cti[index_l0] = 0
 

 plothist, list_energies, xhist_input, yhist_input, bin=hw[0]*2, background=1, color=0, boxplot=0
 plothist, list_energies_cti, xhist, yhist, bin=hw[0]*2, /overplot, color=4, boxplot=0
 
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
 
    outfilenameqdp = 'pixel_cti_'+string_split[0]+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'.qdp'
    openw, lu, outfilenameqdp, /get_lun
    print_array = [transpose(xhist_gt0), transpose(replicate(hw[0],n_elements(xhist_gt0))), transpose(yhist_gt0/(hw[0] * 2 * expo_qdp[0])), transpose(err_cti_kev)]
    printf, lu, print_array
    free_lun, lu

    if keyword_set(key_plot) then begin
       entry_device = !d.name
       plotfilenameps = 'pixel_cti_'+string_split[0]+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color,/tt_font, set_font='Times', font_size=10
       titstring = 'Input [black] + distorted [red] '+ string_split[0]
       plothist, list_energies, xhist, yhist, bin=hw[0]*2, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "Cnts/s/keV", boxplot=0
       plothist, list_energies_cti, bin=hw[0]*2, /overplot, color=3, boxplot=0
       expostr = 'Expo = '+ strtrim(string(floor(expo_qdp[0])), 2)+ ' s'
       xyouts, 0.75, 0.8, expostr, /norm;, charsize=1.2
       driftstr = 'Drift = '+ strtrim(string(floor(dety_max[0]-dety_min[0])), 2)+ ' pixels'
       xyouts, 0.75, 0.75, driftstr, /norm;, charsize=1.2
       trmaxstr = 'DETY max = '+ strtrim(string(floor(dety_max[0])), 2)
       xyouts, 0.75, 0.7, trmaxstr, /norm;, charsize=1.2
              
       device,/close  
       set_plot,entry_device
       
       plotfilenameps = 'pixel_cti_'+string_split[0]+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'_log.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color
       titstring = 'Input [black] + distorted [red] ' + string_split[0]
       plothist, list_energies, bin=hw[0]*2, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "Cnts/s/keV",xr=[0.2,3],/ylog, yr=[1, max(yhist)*1.1], xst=1, yst=1, boxplot=0
       plothist, list_energies_cti, bin=hw[0]*2, /overplot, color=3, boxplot=0
       xyouts, 0.75, 0.8, expostr, /norm;, charsize=1.2
       xyouts, 0.75, 0.75, driftstr, /norm;, charsize=1.2
       xyouts, 0.75, 0.7, trmaxstr, /norm;, charsize=1.2
       device,/close  
       set_plot,entry_device

       ;Additional plots of distribution of losses vs DETY and ENERGY

       !P.Multi = [0, 1, 2, 0, 0]
        
       plotfilenameps = 'pixel_losses_vs_dety_energy_'+string_split[0]+'.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color
       titstring = 'CTI losses '+ string_split[0]
       plot, list_energies, losses_array_ev, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "CTI Losses (eV)",xr=[0.2,3], psym=1, boxplot=0
       plot, n_transfers, losses_array_ev, background=1, color=0, tit = titstring, xtit="DETY", ytit = "CTI Losses (eV)", psym=1, boxplot=0
       
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

 
