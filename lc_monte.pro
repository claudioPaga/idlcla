;@fit_lc
pro lc_monte_wrap,n=n,int=int,nsig=nsig
  
  cd,!mdata
  dir=file_search('GRB*')
  ndir=n_elements(dir)
  g=0
  nw=ndir
  stop
  for i=g,nw-1 do begin
     cd,dir[i]
     print,dir[i],i
     lfile='lc_fit_out_idl_int7.dat'
     if exist(lfile) then begin 
        if numlines(lfile) gt 2 then begin 
;           if exist('flares_gtis.dat') or not exist('lc_fit_out_idl_int7.dat') then $
           lc_monte,/noplot,/ps,nsim=n,int=int,nsig=nsig
        endif 
     endif 
     cd,'..'
  endfor 
  
  return
end 

pro lc_monte,time,tt,exptime,cts,err,corr_fact,srccounts,backcounts,p,mo,pname,outperr,noplot=noplot,ps=ps,nsim=nsim,nsig=nsig,parinfo=parinfo
  
  pt=systime(1)
  if n_elements(nsim) eq 0 then nsim=1000
  print,'n = ',nsim
  intmo='int_'+mo
  
  np=n_elements(p)
  n=nsim

  ;;; 2 sigma confidence level
  if n_elements(nsig) ne 0 then begin 
     if nsig eq 2 then conf=0.9544997 
     if nsig eq 3 then conf=0.997
  endif else conf=0.9
  
  conf1=(1.-conf)/2.
  conf2=(1.-conf)/2.+conf
  print,'Using '+ntostr(conf)+' confidence interval'
  ;;; 90% confidence level
;  conf1=0.05
;  conf2=0.95
  
;  ncounts=exptime*lc.src_rate
  ncounts=srccounts+backcounts
  wb0=where(backcounts eq 0,nwb0)
  if nwb0 gt 0 then backcounts[wb0]=1e-4
  wdet=where(err gt 0,nlc)
;  time=lc.time
  timerr=fltarr(2,nlc)
;  timerr[0,*]=lc.time-lc.tstart
;  timerr[1,*]=lc.tstop-lc.time
  timerr[0,*]=time[wdet]-tt[0,wdet]
  timerr[1,*]=tt[0,wdet]-time[wdet]
  
  tmp=execute('fp='+mo+'(time,p)')
  pp=fltarr(n_elements(p),nsim)
  seed=53464813d
  for i=0,nsim-1 do begin ;; loop over each simulated fit
     if i mod 100 eq 0 then print,'.',format='(a,$)'
     fcounts=fltarr(nlc)
     bcounts=fcounts
     for j=0,nlc-1 do begin ;;loop over each data point
        fcounts[j]=randomn(seed,1.,poi=ncounts[j])
        bcounts[j]=randomn(seed,1.,poi=backcounts[j])
     endfor 
;     frate=fcounts/exptime
;     fraterr=sqrt(fcounts)/exptime
     frate=(fcounts-bcounts)*corr_fact/exptime
     fraterr=sqrt(fcounts+bcounts)*corr_fact/exptime
     
     ;;do fit
     newp=mpfitfun(intmo,tt[*,wdet],frate[wdet],fraterr[wdet],p,yfit=yfit,perror=perror,nprint=100,parinfo=parinfo,dof=dof,/quiet)
     
     if not finite(newp[0]) then newp=p
;     if not keyword_set(noplot) then begin 
;        !p.multi=[0,1,2]
;        ploterror,time[wdet],cts[wdet],timerr[0,wdet],err[wdet],/xlog,/ylog,psym=3,/nohat,xtitle='time',ytitle='counts/s',yrange=yrange
;        oplot,time,fp,color=!green
;     endif 
     
;     tmp=execute('fnewp='+mo+'(t2,newp)')
;     if not keyword_set(noplot) then begin 
;        ploterror,time[wdet],frate[wdet],timerr[0,wdet],fraterr[wdet],/xlog,/ylog,psym=3,/nohat,;xtitle='time',ytitle='counts/s',yrange=yrange
;        oplot,t2,fnewp,color=!green
;        !p.multi=0
;     endif 
     pp[*,i]=newp
;     print,newp,chisq2
;     k=get_kbrd(10)
;     if k eq 's' then stop
  endfor 

;  if keyword_set(ps) then begplot,name='lc_monte_err_plots_int'+int+'.eps',/encap,/land
  nnp=fix(np/2)
  !p.multi=[0,nnp-1,nnp]
;  !x.omargin=[10,0]
;  !x.margin=[8,1]
  outp=fltarr(np)
  outperr=fltarr(2,np)
;  alpha=!tsym.alpha
;  pnames2=['Norm',alpha+'!L1!N','t!Lb,1!N',alpha+'!L2!N','t!Lb,2!N',alpha+'!L3!N','t!Lb,3!N',alpha+'!L4!N']

  for i=0,np-1 do begin
     n5=round(n*conf1)-1
     n95=round(n*conf2)-1
     
     if i mod 2 eq 0 and i ne 0 then begin 
        ppp=alog10(pp[i,*]) 
        bin=0.01
        outp[i]=10^median(ppp)
        sp=ppp[sort(ppp)]
        outperr[0,i]=outp[i]-10^sp[n5]
        outperr[1,i]=10^sp[n95]-outp[i]
     endif else begin 
        ppp=pp[i,*]
        bin=0.01
        outp[i]=median(ppp)
        sp=ppp[sort(ppp)]
        outperr[0,i]=outp[i]-sp[n5]
        outperr[1,i]=sp[n95]-outp[i]
     endelse 
     if i gt 0 then begin 
        w0=where(max(ppp)-min(ppp) gt 0,nw0)
        if nw0 gt 0 then begin 
           plothist,ppp,bin=bin,xtitle=pname[i],charsize=2,ytitle='N'
           oplot,[sp[n5],sp[n5]],[0,n],line=1
           oplot,[sp[n95],sp[n95]],[0,n],line=1
        endif 
     endif 
  endfor 
  !p.multi=0
  if keyword_set(ps) then endplot
  print
  print,p
  print,transpose(perror[0,*])
  print,transpose(perror[1,*])
  print,outp
  print,transpose(outperr[0,*])
  print,transpose(outperr[1,*])
  
;  if not keyword_set(nowrite) then begin 
;     print,'writing out lc_fit_out_idl_int'+int+'.dat'
;     openw,lun,'lc_fit_out_idl_int'+int+'.dat',/get_lun
;     norm=p[0]
;     normerr=outperr[0,*]
;     pnames=pnames[1:*]
;     for i=0,n_elements(p)-2 do begin
;        j=i+1
;        printf,lun,pnames[i]+' '+ntostr(p[j])+' '+ntostr(outperr[0,j])+' '+ntostr(outperr[1,j])
;     endfor
;     printf,lun,'Norm '+ntostr(norm)+' '+ntostr(normerr[0])+' '+ntostr(normerr[1])
;     printf,lun,'Chisq '+ntostr(chisq)
;     printf,lun,'dof '+ntostr(fix(dof))
;     close,lun
;     free_lun,lun
;  endif 
  ptime,systime(1)-pt
  
  return
end 
