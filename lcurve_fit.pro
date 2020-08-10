;infile=The filename (a string) of the grb lightcurve (in fits format)
;start_time=Time from BAT trigger to first XRT frame
;Fits the lightcurve with a power law:  Cnts=a0*Time^(a1)+a2

PRO lcurve_fit, infile, start_time

curve=mrdfits(infile,1,hdr)
time=(curve.x)
cnts=(curve.density)


time_pl=time+start_time                     ;Add the "slew" time from BAT trigger to first XRT frame
a0_guess=cnts[5]/(time_pl[0])^(-1.5)        ;
print,'a0 initial guess= ',a0_guess         ;Makes an initial guess of the a0 parameter using a 

in_param=[a0_guess,-1.5,0]                  ;power law index of -1.5

read,a0_guess,prompt='New a0 estimate? '

;in_param=[10000,-1.5,0]

result=comfit(time_pl,cnts,in_param,/geom,yfit=yfit,sigma=erro)  ;result=[a0,a1,a2]. a1 is the power law index


;SCREEN PLOT

res=strsplit(infile,'.',/extract)
titol='GRB '+res[0]+' lightcurve fit'
set_plot,'X'
plot,time_pl,cnts,/xlog,/ylog,yr=[0.1,20],xr=[10,5000],background=16777215,color=0,xtit='Time',ytit='Cnts/s',tit=titol
bs=strmid(string(result[1]),5,6)
bserr=strmid(string(erro[1]),4,5)
oplot,time_pl,yfit,color=0
xyscritta='Power low index: '+bs+'+/-'+bserr
xyouts,0.2,0.2,xyscritta,/norm,color=0

;END SCREEN PLOT

;POSTSCRIPT PLOT

entry_device=!d.name

set_plot,'PS'
psname=res[0]+'lc_fit.ps'
device,filename=psname
plot,time_pl,cnts,/xlog,/ylog,yr=[0.1,20],xr=[10,5000],background=16777215,color=0,xtit='Time',ytit='Cnts/s',tit=titol
oplot,time_pl,yfit,color=0
xyscritta='Power low index: '+bs+'+/-'+bserr
xyouts,0.1,0.1,xyscritta,/norm,color=0
device,/close_file

set_plot,entry_device

;END OF POSTSCRIPT PLOT


;PRINT PARAMETERS AND FIT VALUES IN AN OUTPUT FILE

print,'Parameters fit: ',result
print,'Parameters error: ',erro

out_write=[transpose(time_pl),transpose(cnts),transpose(yfit)]

openw,nomelogico,'grb_lcurve_corr.txt',/get_lun
printf,nomelogico,'        Time         Counts          Fit'
printf,nomelogico,out_write
printf,nomelogico,'Power law fit data:  Counts=a0*T^(a1)+a2'
printf,nomelogico,'Fit param (a0,a1,a2)= ',result
printf,nomelogico,'Param error= ',erro
free_lun,nomelogico

end
