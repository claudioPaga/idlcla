pro distort_cti_qdp_pixel, qdp_filename, expo_qdp, n_transfers, trap_density_pixel, cti_alpha, cti_beta, key_plot = key_plot, key_serial = key_serial, serial_n_transfers = serial_n_transfers, serial_cti_alpha = serial_cti_alpha, serial_cti_beta = serial_cti_beta
;;;
;;; CP, 28 Nov 2018
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
;;; EXAMPLE:
;;; distort_cti_qdp_pixel, 'ccd2_sxi_puppis_strips_dety3893_expo25200.qdp', 25200., 3893, 1., 0.031, 0.44, key_plot=1, key_serial = 1, serial_n_transfers = 1000, serial_cti_alpha = 0.012, serial_cti_beta = 0.49
;;;  
;;; MODIFICATION HISTORY:
;;;       Written by:     
;;;               Claudio Pagani
;;; Adapted from distort_cti_qdp.pro to include pixel-based treatment
;;; after discussion with David Hall at SMILE calibration workshop   
;;;
;;;
;;; CP, 18/1/2019
;;; V2, Modified to implement serialCTI distortion
;;;

  
 if n_params() lt 6 then begin
    print,"Please enter correct parameters: distort_cti_qdp_pixel, qdp_filename, expo_qdp, n_transfers, trap_density_pixel, cti_alpha, cti_beta"
    print,"Please enter correct parameters: distort_cti_qdp_pixel, 'smile_sxi_casa_10ks_323.qdp', 10000, 100, 1., 0.1, 0.5"
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
  counts_bin = round(counts_kev * hw * 2 * expo_qdp)
  counts_tot = long(total(counts_bin))
  list_energies = dblarr(counts_tot)
  list_energies_cti = dblarr(counts_tot)
  losses_array_ev = dblarr(counts_tot)
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


 ;Pixel-based treatment

 for count_events = 0, n_elements(list_energies)-1 do begin

     ; Apply parallel CTI
    
    eve_energy = list_energies[count_events]
    counter = 0
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

    ;Apply serialCTI (also pixel-based treatment)
    
    if key_serial and eve_energy gt 0 then begin
       serial_counter = 0
       repeat begin
       ;Derive loss in pixel, update event energy
          serial_capture_prob = 1. - exp(-serial_cti_alpha * eve_energy^(1-serial_cti_beta))
          if serial_counter eq 0 then serial_capture_prob_store = [serial_capture_prob] else serial_capture_prob_store = [serial_capture_prob_store, serial_capture_prob]
       ;Sample from binomial distribution
          res = randomn(seed, BINOMIAL = [trap_density_pixel , serial_capture_prob])
          losses_kev = res * 3.65/1000.
          eve_energy -= losses_kev
          serial_counter++
       endrep until (serial_counter eq serial_n_transfers or eve_energy le 0.)
    endif
       
    list_energies_cti[count_events] = eve_energy
 endfor
             
 lost = where(list_energies_cti le 0., xrays_lost_count)
 print, 'Number of events losing all energy: ', xrays_lost_count
 
; Remove negative energies, these events will be lost
 index_positive_energies = where(list_energies_cti gt 0., n_positive)
 if n_positive gt 0 then list_energies_cti = list_energies_cti[index_positive_energies] else begin
    print, 'CTI very high, all events are lost!!!'
    return
 endelse
 
 plothist, list_energies, xhist_input, yhist_input, bin=hw[0]*2, background=1, color=0, box=0
 plothist, list_energies_cti, xhist, yhist, bin=hw[0]*2, /overplot, color=4, box=0
 
 string_split = strsplit(qdp_filename, '.qdp', /regex, /extract)
 if key_serial then outfile_name = string_split[0]+'_pixel_cti_detx'+strtrim(string(serial_n_transfers),2)+'_dety'+strtrim(string(n_transfers),2)+'.qdp' else outfile_name = string_split[0]+'_pixel_cti_dety'+strtrim(string(n_transfers),2)+'.qdp'

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
 
    if key_serial then outfilenameqdp = 'pixel_cti_'+string_split[0]+'_detx'+strtrim(string(serial_n_transfers),2)+'_dety'+strtrim(string(n_transfers),2)+'.qdp' else outfilenameqdp = 'pixel_cti_'+string_split[0]+'_dety'+strtrim(string(n_transfers),2)+'.qdp'
    openw, lu, outfilenameqdp, /get_lun
    printf, lu, [transpose(xhist_gt0), transpose(replicate(hw[0],n_elements(xhist_gt0))), transpose(yhist_gt0/(hw * 2 * expo_qdp)), transpose(err_cti_kev)]
    free_lun, lu

    if keyword_set(key_plot) then begin
       entry_device = !d.name
       if key_serial then plotfilenameps = 'pixel_cti_'+string_split[0]+'_detx'+strtrim(string(serial_n_transfers),2)+'_dety'+strtrim(string(n_transfers),2)+'.ps'  else plotfilenameps = 'pixel_cti_'+string_split[0]+'_dety'+strtrim(string(n_transfers),2)+'.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color,/tt_font, set_font='Times', font_size=10
       if key_serial then plotfilenameps = 'pixel_cti_'+string_split[0]+'_detx'+strtrim(string(serial_n_transfers),2)+'_dety'+strtrim(string(n_transfers),2)+'_log.ps'  else plotfilenameps = 'pixel_cti_'+string_split[0]+'_dety'+strtrim(string(n_transfers),2)+'_log.ps'
       plothist, list_energies, xhist, yhist, bin=hw[0]*2, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "Cnts/s/keV", box = 0, charsize=0.8
       plothist, list_energies_cti, bin=hw[0]*2, /overplot, color=4, box=0
       
       device,/close  
       set_plot,entry_device

       ;Log Version
       
       if key_serial then plotfilenameps = 'pixel_cti_'+string_split[0]+'_detx'+strtrim(string(serial_n_transfers),2)+'_dety'+strtrim(string(n_transfers),2)+'_log.ps'  else plotfilenameps = 'pixel_cti_'+string_split[0]+'_dety'+strtrim(string(n_transfers),2)+'_log.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color
       if key_serial then titstring = 'Input [black] + distorted [blue] - DETX = ' +strtrim(string(serial_n_transfers),2)+ ' DETY = ' +strtrim(string(n_transfers),2)+ ', ' + string_split[0] else titstring = 'Input [black] + distorted [blue] - DETY = ' +strtrim(string(n_transfers),2)+ ', ' + string_split[0]
       plothist, list_energies, bin=hw[0]*2, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "Cnts/s/keV",xr=[0.2,3],/ylog, yr=[1, max(yhist)*1.1], xst=1, yst=1, box = 0, charsize=0.8
       plothist, list_energies_cti, bin=hw[0]*2, /overplot, color=3, box =0
       
       device,/close  
       set_plot,entry_device


       ;Side by side linear version

       !P.Multi = [0, 2, 2, 0, 0]
       entry_device = !d.name
       if key_serial then plotfilenameps = 'pixel_cti_'+string_split[0]+'_detx'+strtrim(string(serial_n_transfers),2)+'_dety'+strtrim(string(n_transfers),2)+'split.ps' else plotfilenameps = 'pixel_cti_'+string_split[0]+'_dety'+strtrim(string(n_transfers),2)+'split.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color,/tt_font, set_font='Times', font_size=10
       plothist, list_energies, xhist, yhist, bin=hw[0]*2, background=1, color=0, xtit="Energy [keV]", ytit = "Input Spectrum (Cnts/s/keV)", box=0, charsize=1.
       plothist, list_energies_cti, bin=hw[0]*2, background=1, color=0, xtit="Energy [keV]", ytit = "CTI Spectrum (Cnts/s/keV)", box=0, charsize=1.
       plothist, list_energies_cti, bin=hw[0]*2, /overplot, color=4, box=0
       
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
stop
; print, 'Losses at 6 keV = ' ,losses6keV
; print, 'Losses at 1 keV = ' ,losses1keV
; print, 'Losses at .5 keV = ',losses500ev
end

 
