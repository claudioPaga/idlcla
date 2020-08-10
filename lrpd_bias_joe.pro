PRO lrpd_bias_joe 

;Reads a catalog file (file.cat) to plot bias values, creates 2 ps
;files, one with the Bias level subtraced on board and the Median Bias
;values calculated from the 20pixels in the header.  The second plot
;is the difference of the two values


readcol,'bias.txt',bias_med,bias_sott,format='(i,i)'
xpp=indgen(115)


index_lrpd=where(bias_med eq 5000)
bias_sott_lrpd=bias_sott(index_lrpd)
xpp_lrpd=xpp(index_lrpd)

index_pupd=where(bias_med ne 5000)
bias_sott_pupd=bias_sott(index_pupd)
xpp_pupd=xpp(index_pupd)

set_plot,'ps'
psname='GRB050117_bias.ps'
device,filename=psname,xs=20,ys=15,/col

tito='GRB050117. Bias subtracted vs Calculated Bias (from 20pix median)'
plot,xpp,bias_med,psym=6,xr=[-5,110],yr=[100,600],xtit='Frame sequence',ytit='BIAS',tit=tito
legenda1='Squares = Bias calculated from median of 20 pixels'
legenda2='!4D!3 = Bias subtraced on board'
xyouts,0.2,0.85,legenda1,/normal
xyouts,0.2,0.8,legenda2,/normal
oplot,xpp_lrpd,bias_sott_lrpd,psym=1
oplot,xpp_pupd,bias_sott_pupd,psym=5

device,/close

end
