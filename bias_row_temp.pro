PRO bias_row_temp, infile, infile2

;Reads a catalog file (file.cat) to plot bias values as a function
;of temperature, creates 2 ps.  The bias is calculated as the median
;of the 20 pixels in the header
;

readcol,infile,frame,format='(a)'

readcol,infile2,frame2,format='(a)'

ccd_t=fltarr(n_elements(frame))
median_bias=fltarr(n_elements(frame))
sigma=fltarr(n_elements(frame))

for i=0,n_elements(frame)-1 do begin
    readframe=frame[i]
    tab=mrdfits(readframe,0,hd0,/silent)
    ccd_t[i]=sxpar(hd0,'T_CCD')
    readframe=frame2[i]
    tab2=mrdfits(readframe,0,hd0,/silent,/dscale)
    median_bias[i]=median(tab2[0:199,*])
    print,tab2[0:199,*]
    sigma[i]=stdev(tab2[0:199,*])
;      if ((ccd_t[i] lt -52) and (bias_med[i] gt 200)) then print,readframe 
endfor

outprint=[transpose(ccd_t(sort(ccd_t))),transpose(median_bias(sort(ccd_t))),transpose(sigma(sort(ccd_t)))]
openw,nomelogico,'ccdtemp_biasrow_median.dat',/get_lun
printf,nomelogico,'CCD temperature and average bias and sigma (from 20 header pixels) for LrPD frames'
printf,nomelogico,'CCD temp       Bias Row Median value       Sigma'
printf,nomelogico,outprint
free_lun,nomelogico

;ave=fltarr(13)
;ccd_ave_plot=intarr(13)
;avesigma=fltarr(13)
;aveindex=0
;for temp_incr=-57,-45 do begin
;    index=where(ccd_t gt temp_incr-0.5 and ccd_t lt temp_incr+0.5)
;    aveindex=(temp_incr+57)
;    ave[aveindex]=mean(bias_med[index])
;    ccd_ave_plot[aveindex]=temp_incr
;    avesigma[aveindex]=mean(sigma[index])
;endfor

!p.multi=[0,1,2,0,0]
set_plot,'ps'
psname=infile2+'_vs_temp.ps'
device,filename=psname,xs=20,ys=15,/col
tito='Bias Row Median value vs CCD temp'
plot,ccd_t,median_bias,psym=5,xtit='CCD temp',ytit='BIAS',tit=tito
plot,ccd_t,sigma,psym=5,xtit='CCD temp',ytit='Sigma'
device,/close
!p.multi=0

end
