nomeps='~/idlplot.ps'
set_plot,'ps'
loadct,39
!P.CHARSIZE = 1.4
!p.thick=7
XYOUTS,.5,.5,'!6'
device,filename=nomeps,xs=15,ys=15,/col        

readcol,'~/tab.txt',x,y1,y2,y3
plot,x,y3,psym=4,sym=2,xtit='Flux',ytit='Photon Index',tit='PILE-UP', $
xr=[0,1010000],yr=[1,3],/xst,/yst
oplot,x,y3,psym=4,sym=2,color=245
oplot,x,y2,psym=4,sym=2,color=145
oplot,x,y1,psym=4,sym=2,color=45
oplot,x,y3,color=245
oplot,x,y2,color=145
oplot,x,y1,color=45



device,/close
spawn,'gv '+nomeps
set_plot,'x'
!p.thick=1
!p.charsize=1
!p.charthick=1
bw


fine:
END
