function const_2gau_lin, X, P

;P=[CONSTANT, CENTROID1, SIGMA1, HEIGHT1, CENTROID2, SIGMA2, HEIGHT2, SLOPE ]

y = dblarr(n_elements(X))

y = P[0] + P[3] * exp(-0.5 * ((X-P[1])/P[2])^2) +  P[6] * exp(-0.5 * ((X-P[4])/P[5])^2)  + P[7] * X

return, y

end

pro fit_2gau_lin_cti_qdp, file_counts_qdp, emin, emax, centroid_guess1, centroid_guess2, key_plot = key_plot

  ;;;CP, 03 Oct 2018
  ;;;
  ;;;Summary - Reads in distorted qdp spectral file (of N132d), fits two
  ;;;          Gaussians + Lin to the brightest emission lines
  ;;;
  ;;;Modified version of fit_2gau_lin_qdp, with explicit identification
  ;;;of the two lines centroids.

  
  if n_params() lt 5 then begin
     print,"Please enter correct parameters: fit_2gau_lin_cti_qdp,,filename,en_min, en_max, centroid_guess1, centroid_guess2,"
     print,"Please enter correct parameters: fit_2gau_lin_cti_qdp,'cti_sxi_n132d_55ks_999.qdp',0.3, 0.45, 0.35, 0.4"
     return
  endif


  DEVICE,DECOMPOSED=0.
  tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1


readcol, file_counts_qdp, en, hw, counts_kev, err_kev, format='(d,d,d,d)'
;Convert the counts/keV to counts in the energy bin (it's
;half-width in the qdp file)
counts_bin = counts_kev * hw * 2
err_bin = err_kev * hw * 2
fit_en_min=emin ;Range of energy for which I create the histogram (that is, the spectrum)
fit_en_max=emax


index_fit = where(en ge fit_en_min and en lt fit_en_max)

en_fit = en[index_fit]
counts_bin_fit = counts_bin[index_fit]
err_bin_fit = err_bin[index_fit]+0.001 ; Added 0.001 of systematic errors
peakcounts=max(counts_bin_fit,peak_index)
plot, en , counts_bin,BACKGROUND=1,color=0, xr=[fit_en_min,fit_en_max],xst=1,psym=1
mincounts=min(counts_bin_fit)
peakenergy = en_fit[peak_index]


;Preliminary fit of the spectrum with a single Gauss, useful mostly to
;determine lines area

;guess_estimates = [total(en_fit*counts_bin_fit), (emax-emin)*0.5, (emax-emin)*0.25, mincounts]
guess_estimates = [peakenergy, (emax+emin)*0.5, (emax-emin)*0.5, mincounts]

yfitgau = GAUSSFIT(en_fit, counts_bin_fit, coeff, MEASURE_ERROR =  err_bin_fit, NTERMS=4, ESTIMATES = guess_estimates)
print,'Coeff of single Gaussian: ', coeff
stop

start = dblarr(8)
;start = [fit_para_counts_bin(3),fit_para_counts_bin(1),4.5,fit_para_counts_bin(0),fit_para_counts_bin(1),4.5,fit_para_counts_bin(0)]

;I'm haveing issues with the fit.
;So I fix and limit some parameters
;I fix the constant (P[0]
;I limit everything else

;Here I will create the array of structures that is used by mpfitfun
;The array has 9 elements (that is, 9 structures), one relative at each
;parameter of the fit
;The structure allows to apply limits to the parameters.
;In particular, it allows a parameter to be fixed or not.
;For example, if we want the second parameter fixed to the value of 666:
;pi(2).fixed=1 (1= fixed, 0=free)
;start(2) = 666
;If we want to limit the third parameter to be within 666 and 777:
;pi(3).limited(0)=1 so the lower limit is limited
;pi(3).limited(1)=1 so the upper limit is limited
;pi(3).limits(0) = 666
;pi(3).limits(1) = 777

pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D], mpprint:0},8)

;I fix the constant to the value derived from the gaussian fit
;pi(0).fixed = 1
start(0) = mincounts

;I limit the centroid
start[1] = centroid_guess1
pi(1).limited(0) = 1
pi(1).limited(1) = 1
pi(1).limits(0) = centroid_guess1-0.1
pi(1).limits(1) = centroid_guess1+0.1

;Limits on the sigmas, the idea is that the column is due to two narrow lines
; left first
start[2] = [0.025]
pi(2).limited(0) = 1
pi(2).limited(1) = 1
pi(2).limits(0) = 0.02
pi(2).limits(1) = 0.035

;I limit the centroid
start[4] = centroid_guess2
pi(4).limited(0) = 1
pi(4).limited(1) = 1
pi(4).limits(0) = centroid_guess2-0.05
pi(4).limits(1) = centroid_guess2+0.05

;Limit the sigma
start[5] = [0.025]
pi(5).limited(0) = 1
pi(5).limited(1) = 1
pi(5).limits(0) = 0.02
pi(5).limits(1) = 0.035


;Limits on lines area (very simple contraints based on Gaussian Fit,
;coeff[0] = Area of Gaussian in fit)
start[3] = coeff[0]
pi(3).limited(0) = 1
pi(3).limited(1) = 1
pi(3).limits(0) = coeff[0]*0.8
pi(3).limits(1) = coeff[0]*1.2
start[6] = coeff[0]
pi(6).limited(0) = 1
pi(6).limited(1) = 1
pi(6).limits(0) = coeff[0]*0.8
pi(6).limits(1) = coeff[0]*1.2

;Limits on the linear splot, it should be positive, given the shape of
;the spectrum at low energies
start[7] = [0.1]

pi(7).limited(0) = 1
;pi(7).limited(1) = 1
pi(7).limits(0) = 0
;pi(7).limits(1) = 1


;I use MPPRINT = 0 so mpfit won't print the parameter values at each iteration 
pi(*).MPPRINT = 0



result = MPFITFUN('const_2gau_lin', en_fit, counts_bin_fit, err_bin_fit, start, PARINFO=pi)
print, result
oplot, en_fit, const_2gau_lin(en_fit, result), linestyle = 3, color = 4, thick=2

;Plot of the single gaussians

g1 = result[0] + result[3] * exp(-0.5 * ((en_fit-result[1])/result[2])^2) + result[7] * en_fit
g2 = result[0] + result[6] * exp(-0.5 * ((en_fit-result[4])/result[5])^2) + result[7] * en_fit
oplot, en_fit, g1, color=3
oplot, en_fit, g2, color=3

;Plot of bb contribution

bb = result[0] + result[7] * en_fit
oplot, en_fit, bb, color=3


;fwhm=2.35*(result(2)+result(3))*0.5*10.  ;With this normalization it's in eV

print,'E_min        Emax     Centroid1(eV)    Sigma1(eV)    Centroid2(eV)   Sigma2(eV)    Delta_centroids'
print,emin*1000.,emax*1000.,result[1]*1000.,result(2)*1000.,result(4)*1000.,result(5)*1000.,(result[4]-result(1))*1000.

string_split = strsplit(file_counts_qdp, '.qdp', /regex, /extract)
outfile_name = string_split[0]+'_qdp_fit_2gau_bb.txt'

openw, lu, outfile_name, /get_lun, WIDTH=250
printf, lu, emin*1000.,emax*1000.,result[1]*1000.,result(2)*1000.,result(4)*1000.,result(5)*1000.,(result[4]-result(1))*1000.
free_lun, lu


if keyword_set(key_plot) then begin
   entry_device = !d.name
   plotfilenameps = string_split[0]+'_qdp_fit_2gau_lin_cti.ps'
   set_plot,'ps'
   device, filename = plotfilenameps, /color,/tt_font, set_font='Times', font_size=10
   titstring = 'N132D fit, CTI distorted spectrum, 2 gau + lin, file = ' + string_split[0]
   plot, en , counts_bin,BACKGROUND=1,color=0, xr=[fit_en_min,fit_en_max],xst=1,psym=1, tit = titstring
   oplot, en_fit, const_2gau_lin(en_fit, result), linestyle = 3, color = 4, thick=2
   oplot, en_fit, g1, color=3
   oplot, en_fit, g2, color=3
   oplot, en_fit, bb, color=3

   device,/close  
   set_plot,entry_device
endif
   
end
