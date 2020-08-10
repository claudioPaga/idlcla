function asym_gau_const, X, P

;P=[BASE, CENTROID, SIGMA1, SIGMA2, HEIGHT]
y = dblarr(n_elements(X))
indexlow = where (X lt P[1])
indexhigh = where (X ge P[1])
y(indexlow) = P[4] * exp(-0.5 * ((X(indexlow)-P[1])/P[2])^2) + P[0]
y(indexhigh) = P[4] * exp(-0.5 * ((X(indexhigh)-P[1])/P[3])^2) + P[0]
return, y

end

pro fit_asym_line,filename,emin,emax,printplotname=printplotname

if n_params() lt 3 then begin 
    print,'syntax - fit_asym_line_lin,evename.evt,enmin,enmax,printplotname=plotname.png'
    return
endif



DEVICE,DECOMPOSED=0.
tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1

tab=mrdfits(filename,1,hd1,/silent)
histobinsize=2
channel_shift=fix(histobinsize*0.5)
fit_en_min=emin ;Range of energy for which I create the histogram (that is, the spectrum)
fit_en_max=emax

colcentral=tab.PI
corrected = HISTOGRAM( colcentral , BINSIZE=histobinsize, MAX=fit_en_max , MIN=fit_en_min)
channels=(indgen(n_elements(corrected)))*histobinsize+fit_en_min
peakcounts=max(corrected,peak_index)
plot,channels,corrected,BACKGROUND=1,color=0,xr=[fit_en_min,fit_en_max],yr=[0,peakcounts+20],tit=title,xst=1,yst=1,psym=1
mincounts=min(corrected)
peakenergy=channels(peak_index)
peak_index_guess=fix(n_elements(corrected)*1./2.)
peak_index_low=fix(peak_index_guess-peak_index_guess*0.9)
peak_index_high=fix(peak_index_guess+peak_index_guess*0.5)
if (peak_index lt peak_index_low) or (peak_index gt peak_index_high) then peak_index=peak_index_guess
fitwidth=fix(10./(histobinsize*1.))
if peak_index-fitwidth gt 0 then begin
    si_en_fit=channels[peak_index-fitwidth:peak_index+fitwidth] 
    si_counts_fit=corrected[peak_index-fitwidth:peak_index+fitwidth]
endif else begin
    si_en_fit=channels[0:peak_index+fitwidth]
    si_counts_fit=corrected[0:peak_index+fitwidth]
endelse
start_guess=[peakcounts*1.,channels(peak_index),5.,mincounts*1.]
;results_corrected=gaussfit(si_en_fit,si_counts_fit,fit_para_corrected,nterms=4,sigma=par_errors,estimates=start_guess,chisq=chivalue)
results_corrected=gaussfit(channels,corrected,fit_para_corrected,nterms=4,sigma=par_errors,estimates=start_guess,chisq=chivalue)
print,fit_para_corrected
oplot,channels,results_corrected,color=4
fwhm_gau = fit_para_corrected[2] * 2.35 * 10. ;FWHM of the gaussian fit in eV

start = [fit_para_corrected(3),fit_para_corrected(1),fit_para_corrected(2),fit_para_corrected(2),fit_para_corrected(0)]

;I'm haveing issues with the fit.
;So I fix and limit some parameters
;I fix the constant (P[0]
;I limit everything else

;Here I will create the array of structures that is used by mpfitfun
;The array has 5 elements (that is, 5 structured), one relative at each
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

pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},5)

;I fix the constant to the value derived from the gaussian fit
pi(0).fixed = 1
start(0) = fit_para_corrected(3)

;I limit the centroid

pi(1).limited(0) = 1
pi(1).limited(1) = 1
pi(1).limits(0) = fit_para_corrected(1)-4
pi(1).limits(1) = fit_para_corrected(1)+4

;Limits on the sigmas, left first

;pi(2).limited(0) = 1
;pi(2).limited(1) = 1
;pi(2).limits(0) = fit_para_corrected(2)-1.
;pi(2).limits(1) = fit_para_corrected(2)+fit_para_corrected(2)*2.
;now to the right
;pi(3).limited(0) = 1
;pi(3).limited(1) = 1
;pi(3).limits(0) = fit_para_corrected(2)-2.
;pi(3).limits(1) = fit_para_corrected(2)+fit_para_corrected(2)*2.

;Also the area should be similar, so I let it vary by 20% max

pi(4).limited(0) = 1
pi(4).limited(1) = 1
pi(4).limits(0) = fit_para_corrected(0)-fit_para_corrected(0)*0.2
pi(4).limits(1) = fit_para_corrected(0)+fit_para_corrected(0)*0.2

;***********END OF PARS LIMIT SETTINGS****************

rerr = sqrt (corrected)
result = MPFITFUN('asym_gau_const', channels, corrected, rerr, start, PARINFO=pi)
print, result

;WARNINGS if best fit parameters hit the limits (there might be something
;wrong in the fit)
if result[1] eq pi(1).limits(0) or result[1] eq pi(1).limits(1) then print,'Centroid hit limit! P(1) =',result[1]
if result[4] eq pi(4).limits(0) or result[4] eq pi(4).limits(1) then print,'Line area hit limit! P(4) =',result[4]


;PLOTS and PRINT results
oplot, channels, asym_gau_const(channels, result), color=3, thick=5
fwhm=2.35*(result(2)+result(3))*0.5*10.  ;With this normalization it's in eV
print,'E_min        Emax     Centroid(eV) Sigmaleft    Sigmaright   FWHM(eV)'
print,emin,emax,result[1]*10.,result(2)*10.,result(3)*10.,fwhm
print,'FWHM of the symmetric gaussian fit =',fwhm_gau

;Hardcopy of the plot if requested

if keyword_set(printplotname) then begin
    nometitolo=' '
    read,nometitolo,prompt='Please type observation name: '
    title=string(nometitolo)+' Line fit - Gaussian and Asymmetric Gaussian fits'
  ;  !P.FONT=1
    ;!X.THICK=1.9
    ;!Y.THICK=1.9
    !P.CHARSIZE=1.2
;    !Y.CHARSIZE=1.4
    ;!P.CHARTHICK=1.9
    plot, channels,corrected,BACKGROUND=1,color=0,xr=[fit_en_min,fit_en_max],tit=title,xst=1,yst=1,psym=3,SYMSIZE=1.5,THICK=2.0,xtit='PI';,yr=[0,peakcounts+20],
    !p.color=0
    oploterror, channels,corrected, rerr, ERRTHICK = 2.0 , psym=3
;            [ /NOHAT, HATLENGTH= , ERRTHICK =, ERRSTYLE=, ERRCOLOR =, 
;              /LOBAR, /HIBAR, NSKIP = , NSUM = , ... OPLOT keywords ]
    oplot, channels,results_corrected,color=4, thick=2
    oplot, channels, asym_gau_const(channels, result), color=3, thick=2
    WRITE_PNG, printplotname , TVRD(TRUE=1), /transparent
endif

end
