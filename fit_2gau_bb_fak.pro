function const_2gau_bb, X, P

;P=[CONSTANT, CENTROID1, SIGMA1, HEIGHT1, CENTROID2, SIGMA2, HEIGHT2, SLOPE ]

y = dblarr(n_elements(X))

y = P[0] + P[3] * exp(-0.5 * ((X-P[1])/P[2])^2) +  P[6] * exp(-0.5 * ((X-P[4])/P[5])^2)  + P[7] * X^2 / (P[8]^4 * [exp(X/P[8])-1])

return, y

end


pro fit_2gau_bb_fak, file_counts_fak, emin, emax, key_plot = key_plot

  ;;;CP, 18 Sept 2018
  ;;;
  ;;;Summary - Reads in qdp spectral file (of E0102), fits two
  ;;;          Gaussians + BB with (A * E^2 / B^4*[exp(E/B)-1]) to the brightest emission lines
  

  
  if n_params() lt 3 then begin
     print,"Please enter correct parameters: fit_2gau_bb_fak,filename,en_min, en_max"
     print,"Please enter correct parameters: fit_2gau_bb_fak,'smile_sxi_casa_10ks_999.fak',0.8, 1.6"
     return
  endif



  DEVICE,DECOMPOSED=0.
  tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1

                                ;Reads in spectrum, transform energies
                                ;in keV, add 0.0025keV to have center
                                ;of bin, SMILE rmf has 5ev channels
  tab = mrdfits(file_counts_fak,1,hd1)
  
  en = tab.CHANNEL * 0.005 + 0.0025 ;
  counts_bin = tab.COUNTS
  err_bin = sqrt(counts_bin)
    
  fit_en_min=emin               ;Range of energy for which I create the histogram (that is, the spectrum)
  fit_en_max=emax

  index_fit = where(en ge fit_en_min and en lt fit_en_max)

  en_fit = en[index_fit]
  counts_bin_fit = counts_bin[index_fit]
  err_bin_fit = err_bin[index_fit]+0.001 ; Added 0.001 of systematic errors
  
  peakcounts=max(counts_bin_fit,peak_index)
  plot, en , counts_bin,BACKGROUND=1,color=0, xr=[fit_en_min,fit_en_max],xst=1,psym=1
  mincounts=min(counts_bin_fit)
  peakenergy = en_fit[peak_index]

  start = dblarr(9)
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D], mpprint:0},9)

;I fix the constant to the value derived from the gaussian fit
;pi(0).fixed = 1
start(0) = mincounts


;I fix the constant to the value derived from the gaussian fit
;pi(0).fixed = 1
start(0) = mincounts

;I limit the centroid
start[1] = 1.05
pi(1).limited(0) = 1
pi(1).limited(1) = 1
pi(1).limits(0) = 1.05-0.05
pi(1).limits(1) = 1.05+0.05

;Limits on the sigmas, the idea is that the column is due to two narrow lines
; left first
start[2] = [0.015]
pi(2).limited(0) = 1
pi(2).limited(1) = 1
pi(2).limits(0) = 0.01
pi(2).limits(1) = 0.03

;I limit the centroid
start[4] = 1.35
pi(4).limited(0) = 1
pi(4).limited(1) = 1
pi(4).limits(0) = 1.35-0.05
pi(4).limits(1) = 1.35+0.05


;Limits on the sigma
start[5] = [0.015]
pi(5).limited(0) = 1
pi(5).limited(1) = 1
pi(5).limits(0) = 0.01
pi(5).limits(1) = 0.03


;BB pars, should be positive
start[7] = [0.1]
start[8] = [0.2]

pi(7).limited(0) = 1
;pi(7).limited(1) = 1
pi(7).limits(0) = 0
;pi(7).limits(1) = 1

pi(8).limited(0) = 1
;pi(8).limited(1) = 1
pi(8).limits(0) = 0.
;pi(8).limits(1) = 2.

;I use MPPRINT = 0 so mpfit won't print the parameter values at each iteration 
pi(*).MPPRINT = 0

result = MPFITFUN('const_2gau_bb', en_fit, counts_bin_fit, err_bin_fit, start, PARINFO=pi)
print, result
oplot, en_fit, const_2gau_bb(en_fit, result), linestyle = 3, color = 4, thick=2

;Plot of the single gaussians

g1 = result[0] + result[3] * exp(-0.5 * ((en_fit-result[1])/result[2])^2) + result[7] * en_fit^2 / (result[8]^4 * (exp(en_fit/result[8])-1))
g2 = result[0] + result[6] * exp(-0.5 * ((en_fit-result[4])/result[5])^2) + result[7] * en_fit^2 / (result[8]^4 * (exp(en_fit/result[8])-1))
oplot, en_fit, g1, color=3
oplot, en_fit, g2, color=3

;Plot of bb contribution

bb = result[0] + result[7] * en_fit^2 / (result[8]^4 * (exp(en_fit/result[8])-1))
oplot, en_fit, bb, color=4

;fwhm=2.35*(result(2)+result(3))*0.5*10.  ;With this normalization it's in eV

print,'E_min        Emax     Centroid1(eV)    Sigma1(eV)    Centroid2(eV)   Sigma2(eV)    Delta_centroids'
print,emin*1000.,emax*1000.,result[1]*1000.,result(2)*1000.,result(4)*1000.,result(5)*1000.,(result[4]-result(1))*1000.

string_split = strsplit(file_counts_fak, '.fak', /regex, /extract)
outfile_name = string_split[0]+'_fak_fit_2gau_bb.txt'

openw, lu, outfile_name, /get_lun, WIDTH=250
printf, lu, emin*1000.,emax*1000.,result[1]*1000.,result(2)*1000.,result(4)*1000.,result(5)*1000.,(result[4]-result(1))*1000.
free_lun, lu

if keyword_set(key_plot) then begin
   entry_device = !d.name
   plotfilenameps = string_split[0]+'_fak_fit_2gau_bb.ps'
   set_plot,'ps'
   device, filename = plotfilenameps, /color,/tt_font, set_font='Times', font_size=10
   titstring = 'Cas A fit, 2 gau + bb, file = ' + string_split[0]
   plot, en , counts_bin,BACKGROUND=1,color=0, xr=[fit_en_min,fit_en_max],xst=1,psym=1, tit = titstring
   oplot, en_fit, const_2gau_bb(en_fit, result), linestyle = 3, color = 4, thick=2

   device,/close  
   set_plot,entry_device
endif
   
end
