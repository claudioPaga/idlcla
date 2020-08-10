PRO lrpd_bias_plot, infile 

;Reads a catalog file (file.cat) fo fits LrPD frames, creates 2 ps
;files, one with the Bias level subtraced on board and the Median Bias
;values calculated from the 20pixels in the header.  The second plot
;is the difference of the two values

readcol,infile,frame,format='(a)'
forprint, frame

bias_sub=intarr(n_elements(frame))
bias_real=intarr(n_elements(frame))
bias_diff=intarr(n_elements(frame))
frame_num=intarr(n_elements(frame))

for i=0,n_elements(frame)-1 do begin
    readframe=frame[i]
    tab=mrdfits(readframe,0,hd0,/silent)
    frame_num[i]=sxpar(hd0,'CCD_FRAM')
    bias_sub[i]=sxpar(hd0,'BIAS_LVL')
    bias_real[i]=sxpar(hd0,'MEDIANPX')
    bias_diff[i]=bias_real[i]-bias_sub[i]
endfor

outprint=[transpose(frame_num),transpose(bias_sub),transpose(bias_real),transpose(bias_diff)]
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
oplot,frame_num,bias_real,psym=6

device,/close

tito=strmid(infile,0,17)+' Bias difference'
set_plot,'ps'
psname=infile+'_bias_diff.ps'
device,filename=psname,xs=15,ys=15,/col
plot,frame_num,bias_diff,psym=6,xtit='FRAME NUMBER',ytit='BIAS DIFFERENCE',tit=tito
device,/close

end
