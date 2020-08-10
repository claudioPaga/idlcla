PRO lrpd_bias_corr, infile 

;Reads a catalog file of fits LrPD frames and subtract the correct
;bias calculated as the difference of the Median value of the 
;20 pixels in the header and the bias level already subtracted
;Creates fits files with updated PHA values with the new bias calculation

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
    print,bias_diff[i]
    tt=mrdfits(readframe,1,hdr,/silent)
    tt.pha=tt.pha-bias_diff[i]
    writefile=frame[i]+'_bias_corr.fits'
    mwrfits,tab,writefile,hd0
    mwrfits,tt,writefile,hdr
endfor
print,'Bias subrtacted'
print,bias_sub
print,'Bias median'
print,bias_real

end
