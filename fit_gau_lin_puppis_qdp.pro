function const_gau_lin, X, P

;P=[CONSTANT, CENTROID1, SIGMA1, HEIGHT1, SLOPE ]

y = dblarr(n_elements(X))

y = P[0] + P[3] * exp(-0.5 * ((X-P[1])/P[2])^2) + P[4] * X

return, y

end

pro fit_gau_lin_puppis_qdp, file_counts_qdp, emin, emax, centroid_guess1, key_plot = key_plot

  ;;;CP, 22 Oct 2018
  ;;;
  ;;;Summary - Reads in qdp spectral file (of Puppis A), fits
  ;;;          Gaussians + Lin to the Silicon emission lines
  ;;;
  ;;;Modified version of fit_2gau_lin_qdp, with explicit identification
  ;;;of the line centroid.
  ;;;
  ;;;Modified version of fit_2gau_lin_cti_qdp, to try to improve line fitting
  

  
  if n_params() lt 4 then begin
     print,"Please enter correct parameters: fit_gau_lin_puppis_qdp,filename,en_min, en_max, centroid_guess1, key_plot = 1"
     print,"Please enter correct parameters: fit_gau_lin_puppis__qdp,'sxi_puppis_30ks_999.qdp',1.7, 2.0, 1.85"
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
guess_estimates = [peakenergy, (emax+emin)*0.5, (emax-emin)*0.5, 0., 4.]

;yfitgau = GAUSSFIT(en_fit, counts_bin_fit, coeff, MEASURE_ERROR =  err_bin_fit, NTERMS=5, ESTIMATES = guess_estimates)
;print,'Coeff of single Gaussian: ', coeff
;stop

start = dblarr(5)
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

pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D], mpprint:0},5)

;I fix the constant to the value derived from the gaussian fit
;pi(0).fixed = 1
start(0) = counts_bin_fit[0]

;I limit the centroid
start[1] = centroid_guess1
pi(1).limited(0) = 1
pi(1).limited(1) = 1
pi(1).limits(0) = centroid_guess1-0.2
pi(1).limits(1) = centroid_guess1+0.2

;Limits on the sigmas, the idea is that the column is due to two narrow lines
; left first
start[2] = [0.025]
pi(2).limited(0) = 1
pi(2).limited(1) = 1
pi(2).limits(0) = 0.015
pi(2).limits(1) = 0.055

;Limits on lines normal factor
;Area_lines = Total_counts - Estimate of area under linear component
;Gauss_factor = Area_lines/(sigma * sqrt(2*3.14))

area_linear = 0.5 * (counts_bin_fit[n_elements(counts_bin_fit)-1] + counts_bin_fit[0]) * (en_fit[n_elements(en_fit)-1]-en_fit[0])
cts_tot = total(counts_bin_fit) * (en_fit[1]-en_fit[0])
area_lines = cts_tot - area_linear
normal_guess = area_lines/(sqrt(2.*!PI)*start[2])
start[3] = normal_guess
pi(3).limited(0) = 1
pi(3).limited(1) = 1
pi(3).limits(0) = start[3]*0.3
pi(3).limits(1) = start[3]*2.

;Limits on the linear slope, it should be negative, given the shape of
;the spectrum around Si
start[4] = [-0.5]
;pi(4).limited(0) = 1
pi(4).limited(1) = 1
pi(4).limits(1) = 0
;pi(7).limits(1) = 1


;I use MPPRINT = 0 so mpfit won't print the parameter values at each iteration 
pi(*).MPPRINT = 0


result = MPFITFUN('const_gau_lin', en_fit, counts_bin_fit, err_bin_fit, start, PARINFO=pi)
print, 'Initial guesses:'
print, '       Constant         En1              Sigma1          Norm1            Slope'
print, start
print, 'Best fit params:'
print, result
oplot, en_fit, const_gau_lin(en_fit, result), linestyle = 3, color = 4, thick=2


;Warns if parameters have hit constraint
indexlow = where(result le pi.limits(0) and pi.limited(0) ne 0, nhitslow)
if nhitslow gt 0 then print, 'Warning, lower limit hit on parameter(s) ', indexlow, result[indexlow], pi[indexlow].limits[0]

indexhigh = where(result ge pi.limits(1) and pi.limited(1) ne 0, nhitshigh)
if nhitshigh gt 0 then print, 'Warning, upper limit hit on parameter(s) ', indexhigh, result[indexhigh], pi[indexhigh].limits[1]

;Plot of the gaussian

g1 = result[0] + result[3] * exp(-0.5 * ((en_fit-result[1])/result[2])^2) + result[4] * en_fit
oplot, en_fit, g1, color=3

;Plot of linear contribution

linear = result[0] + result[4] * en_fit
oplot, en_fit, linear, color=3


;fwhm=2.35*(result(2)+result(3))*0.5*10.  ;With this normalization it's in eV

print,'        E_min        Emax     Centroid1(eV)    Sigma1(eV)'
print,emin*1000.,emax*1000.,result[1]*1000.,result(2)*1000.

string_split = strsplit(file_counts_qdp, '.qdp', /regex, /extract)
outfile_name = string_split[0]+'_qdp_fit_gau_lin_cti.txt'

openw, lu, outfile_name, /get_lun, WIDTH=250
printf, lu, emin*1000.,emax*1000.,result[1]*1000.,result(2)*1000.
free_lun, lu


if keyword_set(key_plot) then begin
   entry_device = !d.name
   plotfilenameps = string_split[0]+'_qdp_fit_gau_lin.ps'
   set_plot,'ps'
   device, filename = plotfilenameps, /color,/tt_font, set_font='Times', font_size=10
   titstring = 'Puppis A fit, observed spectrum, no CTI, gau + lin, file = ' + string_split[0]
   plot, en , counts_bin,BACKGROUND=1,color=0, xr=[fit_en_min,fit_en_max],xst=1,psym=1, tit = titstring
   oplot, en_fit, const_gau_lin(en_fit, result), linestyle = 3, color = 4, thick=2
   oplot, en_fit, g1, color=3
   oplot, en_fit, linear, color=3

   device,/close  
   set_plot,entry_device
endif

end
