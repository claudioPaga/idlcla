PRO bias_map_temp, infile

;Reads a catalog file (file.cat) to plot bias values as a function
;of temperature, creates 2 ps.  The bias is calculated as the median
;of the 20 pixels in the header
;

readcol,infile,frame,format='(a)'

ccd_t=fltarr(n_elements(frame))
median_bias=fltarr(n_elements(frame))
sigma=fltarr(n_elements(frame))

for i=0,n_elements(frame)-1 do begin
    readframe=frame[i]
    tab=mrdfits(readframe,0,hd0,/silent)
    ccd_t[i]=sxpar(hd0,'T_CCD')
    tab2=tab[180:280,200:400]
    median_bias[i]=median(tab2)
    sigma[i]=stdev(tab2)
endfor

outprint=[transpose(ccd_t(sort(ccd_t))),transpose(median_bias(sort(ccd_t))),transpose(sigma(sort(ccd_t)))]
openw,nomelogico,'ccdtemp_biasmap_median.dat',/get_lun
printf,nomelogico,'Bias Map: CCD temp, median and sigma'
printf,nomelogico,'CCD temp       Bias Map Median       Sigma'
printf,nomelogico,outprint
free_lun,nomelogico

!p.multi=[0,1,2,0,0]
set_plot,'ps'
psname=infile+'_vs_temp.ps'
device,filename=psname,xs=20,ys=15,/col
tito='Bias Map Median value vs CCD temp'
plot,ccd_t,median_bias,psym=5,xtit='CCD temp',ytit='BIAS',tit=tito
plot,ccd_t,sigma,psym=5,xtit='CCD temp',ytit='Sigma'
device,/close
!p.multi=0

end
