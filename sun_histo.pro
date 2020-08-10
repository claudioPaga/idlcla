PRO sun_histo, sunfile, slew_flag_file, plotname 
;
;sunfile:  sun angles file from DTAS
;slew_flag_file:  SAC_MODE channel in DTAS
;plotname:  plot.ps output file
;
;Selects times when SETTLED
;Corrects for 1/sin(Sun angle)  !!!!!! change to program needed
;Plot the INTEGRATED sun angle fraction time distibution,
;the fraction of time spent at angles > of that Sun angle

readcol,sunfile,tsun,sunangle,format='(a,i)'
readcol,slew_flag_file,tsl,value,format='(a,i)'

met_tsl=strarr(n_elements(tsl))
met_tsun=strarr(n_elements(tsun))

for iii=0l,n_elements(tsl)-1 do met_tsl[iii]=date2met(tsl[iii])
for ii=0l,n_elements(tsun)-1 do met_tsun[ii]=date2met(tsun[ii])

sun_inter = INTERPOL(sunangle,double(met_tsun),double(met_tsl),/SPLINE)

index_set=where(value eq 4)
sun_plot_set=sun_inter[index_set]

nomeps=plotname
set_plot,'PS'
device,filename=nomeps

ttt=sun_plot_set
max=185
min=45
binsize=5
stra='SUN ANGLE'
strc='SUN ANGLE INTEGRATED DISTRIBUTION DAY 060-100 - SETTLED'
    nn=fix((max-min)/binsize)
    xas=min+(findgen(nn))*binsize+binsize/2.
    hh=histogram(ttt,min=min,max=max,binsize=binsize)
    hma=1.1*max(hh)
    print,xas
    print,hh
    yplot=hh;/sin(((xas+binsize/2-0.5)/360)*6.28)         ;!!!! if sin correction wanted, modify procedure here!!!
    yplot_sum=total(yplot)
    yplot_frac=yplot/yplot_sum
    ntot=n_elements(yplot)-1
    yplot_int=fltarr(ntot+1)
    for pl_cnt=ntot,0,-1 do yplot_int[pl_cnt]=total(yplot_frac[pl_cnt:ntot])
    print,yplot_int
    plot,xas,yplot_int,psym=10,xtit=stra,tit=strc,yr=[0,1]

device,/close
spawn,'gv '+nomeps

set_plot,'x'
end
