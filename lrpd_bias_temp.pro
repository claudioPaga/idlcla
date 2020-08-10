PRO lrpd_bias_temp, infile

;Reads a catalog file (file.cat) to plot bias values as a function
;of temperature, creates 2 ps.  The bias is calculated as the median
;of the 20 pixels in the header
;

readcol,infile,frame,format='(a)'
;forprint, frame

timestart=fltarr(n_elements(frame))
saa_flag=intarr(n_elements(frame))
bias_sub=intarr(n_elements(frame))
bias_real=intarr(n_elements(frame))
bias_diff=intarr(n_elements(frame))
bias_med=intarr(n_elements(frame))
frame_num=intarr(n_elements(frame))
ccd_t=fltarr(n_elements(frame))
sigma=fltarr(n_elements(frame))
bias_pixel=intarr(20)

for i=0,n_elements(frame)-1 do begin
    readframe=frame[i]
    tab=mrdfits(readframe,0,hd0,/silent)
    saa_flag[i]=sxpar(hd0,'SAA')
    ccd_t[i]=sxpar(hd0,'T_CCD')
    bias_sub[i]=sxpar(hd0,'BIAS_LVL')
    for ii=0,19 do begin
        pix_num=string(ii)
        pixel_name='PIXEL'+strtrim(pix_num,2)
        bias_pixel[ii]=sxpar(hd0,pixel_name)
    endfor
    index_median=where(bias_pixel lt 400,n_pix_nosaa)
    if (n_pix_nosaa ne 0) then begin
        bias_median_values=bias_pixel(index_median)
        bias_med[i]=median(bias_median_values)
        sigma[i]=stdev(bias_median_values)
    endif else bias_med[i]=0
;      if ((ccd_t[i] lt -52) and (bias_med[i] gt 200)) then print,readframe 
endfor

outprint=[transpose(ccd_t),transpose(bias_med)]
openw,nomelogico,'ccdtemp_bias.dat',/get_lun
printf,nomelogico,'CCD temp, Bias,'
printf,nomelogico,outprint
free_lun,nomelogico

res=poly_fit(ccd_t,bias_med,2)
print,'Fit of BIAS=a+b*CCD_T+c*CCD_t^2'
print,'a=',res[0]
print,'b=',res[1]
print,'c=',res[2]

ave=fltarr(13)
ccd_ave_plot=intarr(13)
avesigma=fltarr(13)
aveindex=0
for temp_incr=-57,-45 do begin
    index=where(ccd_t gt temp_incr-0.5 and ccd_t lt temp_incr+0.5)
    aveindex=(temp_incr+57)
    ave[aveindex]=mean(bias_med[index])
    ccd_ave_plot[aveindex]=temp_incr
    avesigma[aveindex]=mean(sigma[index])
endfor

set_plot,'ps'
psname=infile+'_bias_onboard.ps'
device,filename=psname,xs=20,ys=15,/col
tito='Bias level subtracted on board vs CCD temp'
plot,ccd_t,bias_sub,psym=5,xtit='CCD temp',ytit='BIAS',tit=tito
device,/close

!p.multi=[0,1,4,0,0]
set_plot,'ps'
psname=infile+'_bias_calc.ps'
device,filename=psname,xs=20,ys=15,/col
tito='Bias level as median of 20 pixel (without thresh) vs CCD temp'
plot,ccd_t,bias_med,psym=5,xtit='CCD temp',ytit='BIAS',tit=tito
plot,ccd_t,sigma,xtit='CCD temp',ytit='Sigma',psym=4
plot,ccd_ave_plot,ave,xtit='CCD temp',ytit='Ave Bias',psym=6
plot,ccd_ave_plot,avesigma,xtit='CCD temp',ytit='Ave Sigma',psym=6
device,/close
!p.multi=0

openw,nomelogico,'ccdtemp_bias.dat',/get_lun
tabprint=[transpose(ccd_ave_plot),transpose(ave),transpose(avesigma)]
printf,nomelogico,'CCD temperature and average bias and sigma (from 20 header pixels) for LrPD frames'
printf,nomelogico,'CCD temp          Bias           Sigma'
printf,nomelogico,tabprint
free_lun,nomelogico
end
