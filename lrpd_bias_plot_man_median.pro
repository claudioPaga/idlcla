PRO lrpd_bias_plot_man_median, infile 

;Reads a catalog file (file.cat) fo fits LrPD frames, creates 2 ps
;files, one with the Bias level subtraced on board and the Median Bias
;values calculated from the 20pixels in the header.  The second plot
;is the difference of the two values

readcol,infile,frame,format='(a)'
forprint, frame

bias_sub=intarr(n_elements(frame))
bias_real=intarr(n_elements(frame))
bias_diff=intarr(n_elements(frame))
bias_med=intarr(n_elements(frame))
frame_num=intarr(n_elements(frame))
bias_pixel=intarr(20)

for i=0,n_elements(frame)-1 do begin
    readframe=frame[i]
    tab=mrdfits(readframe,0,hd0,/silent)
    frame_num[i]=sxpar(hd0,'CCD_FRAM')
    bias_sub[i]=sxpar(hd0,'BIAS_LVL')
    bias_pixel[0]=sxpar(hd0,'PIXEL0')
    bias_pixel[1]=sxpar(hd0,'PIXEL1')
    bias_pixel[2]=sxpar(hd0,'PIXEL2')
    bias_pixel[3]=sxpar(hd0,'PIXEL3')
    bias_pixel[4]=sxpar(hd0,'PIXEL4')
    bias_pixel[5]=sxpar(hd0,'PIXEL5')
    bias_pixel[6]=sxpar(hd0,'PIXEL6')
    bias_pixel[7]=sxpar(hd0,'PIXEL7')
    bias_pixel[8]=sxpar(hd0,'PIXEL8')
    bias_pixel[9]=sxpar(hd0,'PIXEL9')
    bias_pixel[10]=sxpar(hd0,'PIXEL10')
    bias_pixel[11]=sxpar(hd0,'PIXEL11')
    bias_pixel[12]=sxpar(hd0,'PIXEL12')
    bias_pixel[13]=sxpar(hd0,'PIXEL13')
    bias_pixel[14]=sxpar(hd0,'PIXEL14')
    bias_pixel[15]=sxpar(hd0,'PIXEL15')
    bias_pixel[16]=sxpar(hd0,'PIXEL16')
    bias_pixel[17]=sxpar(hd0,'PIXEL17')
    bias_pixel[18]=sxpar(hd0,'PIXEL18')
    bias_pixel[19]=sxpar(hd0,'PIXEL19')
    index_median=where(bias_pixel lt 2000)
    bias_median_values=bias_pixel(index_median)
    bias_med[i]=median(bias_median_values)
    
;    bias_real[i]=sxpar(hd0,'MEDIANPX')
;    bias_diff[i]=bias_real[i]-bias_sub[i]
     bias_diff[i]=bias_med[i]-bias_sub[i]
endfor

outprint=[transpose(frame_num),transpose(bias_sub),transpose(bias_med),transpose(bias_diff)]
print,'Frame number, Subtracted Bias, Bias Median, Bias difference '
print,outprint

set_plot,'ps'
psname=infile+'_bias_comp.ps'
device,filename=psname,xs=15,ys=15,/col

tito=strmid(infile,0,17)+' Bias subtracted vs Calculated Bias (from 20pix median)'
plot,frame_num,bias_sub,psym=5,yr=[100,400],xtit='FRAME NUMBER',ytit='BIAS',tit=tito
legenda1='Squares = Bias calculated from median of 20 pixels'
legenda2='!4D!3 = Bias subtraced on board'
xyouts,0.2,0.85,legenda1,/normal
xyouts,0.2,0.8,legenda2,/normal
oplot,frame_num,bias_med,psym=6

device,/close

tito=strmid(infile,0,17)+' Bias difference'
set_plot,'ps'
psname=infile+'_bias_diff.ps'
device,filename=psname,xs=15,ys=15,/col
plot,frame_num,bias_diff,psym=6,xtit='FRAME NUMBER',ytit='BIAS DIFFERENCE',tit=tito
device,/close

end
