PRO lrpd_bias_plot_man, infile 

;Reads a catalog file (file.cat) fo fits LrPD frames, creates 2 ps
;files, one with the Bias level subtraced on board and the Median Bias
;values calculated from the 20pixels in the header.  The second plot
;is the difference of the two values

readcol,infile,frame,format='(a)'
;forprint, frame

saa_flag=intarr(n_elements(frame))
bias_sub=intarr(n_elements(frame))
bias_real=intarr(n_elements(frame))
bias_diff=intarr(n_elements(frame))
bias_med=intarr(n_elements(frame))
frame_num=intarr(n_elements(frame))
bias_pixel=intarr(20)

for i=0,n_elements(frame)-1 do begin
    readframe=frame[i]
    tab=mrdfits(readframe,0,hd0,/silent)
    saa_flag[i]=sxpar(hd0,'SAA')
    frame_num[i]=sxpar(hd0,'CCD_FRAM')
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
    endif else bias_med[i]=0
;    bias_real[i]=sxpar(hd0,'MEDIANPX')
;    bias_diff[i]=bias_real[i]-bias_sub[i]
    bias_diff[i]=bias_med[i]-bias_sub[i]
endfor

outprint=[transpose(frame_num),transpose(bias_sub),transpose(bias_med),transpose(bias_diff),transpose(saa_flag)]
print,'Frame number, Subtracted Bias, Bias Median, Bias difference, SAA_Flag '
print,outprint

xplot=indgen(n_elements(frame))



saaplot_index=where(saa_flag eq 1, counts)
saaplot_y=bias_med(saaplot_index)
saaplot_x=xplot(saaplot_index)
print,saaplot_x
print,saaplot_y



set_plot,'ps'
psname=infile+'_bias_comp.ps'
device,filename=psname,xs=20,ys=15,/col

tito=strmid(infile,0,10)+'Bias subtracted vs Calculated Bias (from 20pix median)'
plot,xplot,bias_sub,psym=5,xr=[-1,max(xplot)+1],yr=[-1,600],xtit='Frame sequence',ytit='BIAS',tit=tito
legenda1='Squares = Bias calculated from median of 20 pixels'
legenda2='!4D!3 = Bias subtraced on board'
xyouts,0.2,0.85,legenda1,/normal
xyouts,0.2,0.8,legenda2,/normal
oplot,xplot,bias_med,psym=6
oplot,saaplot_x,saaplot_y,psym=7

device,/close

tito=strmid(infile,0,17)+' Bias difference'
set_plot,'ps'
psname=infile+'_bias_diff.ps'
device,filename=psname,xs=20,ys=15,/col
plot,xplot,bias_diff,psym=6,xtit='Frame sequence',ytit='BIAS DIFFERENCE',tit=tito
device,/close

end
