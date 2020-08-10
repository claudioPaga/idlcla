@fit_functions
pro fit_pow_model,t,cts,terr,err,p,model,pnames,yfit,newp,perror,chisq,dof,weights,src,back,status=status,name=name,silent=silent
  
  if keyword_set(silent) then quiet=1 else quiet=0
  np=n_elements(p)
  parinfo = parinfo_struct(np)

  parinfo.parname=pnames
  parinfo.value=p
  parinfo[0].limits=[0,0] ;;norm > 0
  parinfo[0].limited=[1,0]
  mint=min(t)*1.1
  maxt=max(t)*0.9
;;pow 1 limits > 0
  parinfo[1].limits=[-2,0]
  parinfo[1].limited=[1,0]
  pmargin=2.
  
  case np of 
     4: begin                   ;bknpow
        ;;break limits
        parinfo[2].limits=[mint,maxt]
        parinfo[2].limited=[1,1]
        parinfo[2].mpminstep=1.
        ;;pow limits
        parinfo[3].limits=[p[3]-pmargin,p[3]+pmargin]
        parinfo[3].limited=[1,1]
     end
     6: begin                   ;bkn2pow
        ;;break limits
        parinfo[2].limits=[mint,p[4]]
        parinfo[4].limits=[p[2],maxt]
        parinfo[[2,4]].limited=[1,1]
        parinfo[[2,4]].mpminstep=1.
        ;;pow limits
        parinfo[3].limits=[p[3]-pmargin,p[3]+pmargin]
        parinfo[5].limits=[p[5]-pmargin,p[5]+pmargin]
        parinfo[[3,5]].limited=[1,1]
     end
     8: begin                   ;bkn3pow
        ;;break limits
        parinfo[2].limits=[mint,p[4]]
        parinfo[4].limits=[p[2],p[6]]
        parinfo[6].limits=[p[4],maxt]
        parinfo[[2,4,6]].limited=[1,1]
        parinfo[[2,4,6]].mpminstep=1.
        ;;pow limits
        parinfo[3].limits=[p[3]-pmargin,p[3]+pmargin]
        parinfo[5].limits=[p[5]-pmargin,p[5]+pmargin]
        parinfo[7].limits=[p[7]-pmargin,p[7]+pmargin]
        parinfo[[3,5,7]].limited=[1,1]
     end
     10: begin                  ;bkn4pow
        ;;break limits
        parinfo[2].limits=[mint,p[4]]
        parinfo[4].limits=[p[2],p[6]]
        parinfo[6].limits=[p[4],p[8]]
        parinfo[8].limits=[p[6],maxt]
        parinfo[[2,4,6,8]].limited=[1,1]
        parinfo[[2,4,6,8]].mpminstep=1.
        ;;pow limits
        parinfo[3].limits=[p[3]-pmargin,p[3]+pmargin]
        parinfo[5].limits=[p[5]-pmargin,p[5]+pmargin]
        parinfo[7].limits=[p[7]-pmargin,p[7]+pmargin]
        parinfo[9].limits=[p[9]-pmargin,p[9]+pmargin]
        parinfo[[3,5,7,9]].limited=[1,1]
     end
     else: begin
     end 
  endcase 

;comberr=sqrt((terr/t)^2+(err/cts)^2)
  comberr=err
  weights=1.                    ;/err^2;comberr
  
  tt=dblarr(2,n_elements(t))
  if not keyword_set(phil) then begin 
     tt[0,*]=t-terr
     tt[1,*]=t+terr
  endif else begin
     tt[0,*]=t-terr[0,*]
     tt[1,*]=t+terr[1,*]
  endelse 
  newp=mpfitfun(model,tt,cts,err,p,parinfo=parinfo,$
                bestnorm=chisq,dof=dof,niter=niter,errmsg=errmsg,$
                perror=perror,yfit=yfit,status=status,nprint=10,$
                ftol=1e-15,xtol=1e-15,gtol=1e-25,quiet=quiet)
  chisq=total(((yfit-cts)/err)^2)

  if not keyword_set(silent) then begin
     print,status
     case status of
        0: print,'Improper input'
        1: print,'both actual and predicted relative reduction in sums of squares < FTOL'
        2: print,'relative error between 2 consecutive iterates is < XTOL'
        3: print,'conditions for STATUS = 1 & 2 both hold'
        4: print,'Cosine of angle between fvec and any column of J is || < GTOL'
        5: print,'Maximum iterations'
        6: print,'FTOL too small, no further reduction possible'
        else: print
     endcase
  endif 
  if status ne 0 then $
     perror=perror*sqrt(chisq/dof)
  
  
;weights=1./comberr^2
;chisq=total((cts-yfit)^2*abs(1./comberr^2))

  return
end 
pro plot_lcfit_results,lc,newp,perror,chisq,dof,breaks,leg,pnames,charsize=charsize,noleg=noleg,noerr=noerr,name=name,ps=ps,nocolor=nocolor,xtitle=xtitle,nolines=nolines,_extra=_extra,flux=flux
  
  orange=!p.color & green=orange
  if not keyword_set(nocolor) then begin
     orange=!orange
     if keyword_set(ps) then green=!black else green=!green
  endif   
  if not keyword_set(noerr) then noerr=0
  if keyword_set(ps) then sym=!tsym else sym=!vsym
  time=lc.time
  cts=lc.src_rate*flux
  err=lc.src_rate_err*flux
  tstart=lc.tstart
  tstop=lc.tstop
  t=dblarr(2,n_elements(time)+breaks)
  
  norm=newp[0]
  normerr=perror[*,0]
  pow1=newp[1]
  pow1err=perror[*,1]
  case breaks of
     0: begin 
        t=time
;        t[0,*]=tstart
;        t[1,*]=tstop
        y=pow(t,newp)
        pnames='Pow1'
        leg='Pow1'+' = '+sigfig(pow1,3)
        if not noerr then leg=leg+' !S!E+'+sigfig(pow1err[1],3)+' !R!I-'+sigfig(pow1err[0],3)
     end
     1: begin 
        pow2=newp[3]
        pow2err=perror[*,3]
        break1=newp[2]
        break1err=perror[*,2]
        w1=where(time gt 0 and time lt break1,nw1)
        w2=where(time ge break1)
;        t[0,*]=[tstart[w1],break1,tstart[w2]]
;        t[1,*]=[tstop[w1],break1,tstop[w2]]
        t=[time[w1],break1,time[w2]]
        y=bknpow(t,newp)
        pnames=['Pow1','Breaktime','Pow2']
        leg=['Pow1'+' = '+sigfig(pow1,3), $
             'Breaktime = '+sigfig(break1,5),$
             'Pow2'+' = '+sigfig(pow2,3)]
        
        if not noerr then $
           leg=leg+[' !S!E+'+sigfig(pow1err[1],3)+' !R!I-'+sigfig(pow1err[0],3),$
                    ' !S!E+'+sigfig(break1err[1],5)+' !R!I-'+sigfig(break1err[0],5),$
                    ' !S!E+'+sigfig(pow2err[1],3)+' !R!I-'+sigfig(pow2err[0],3)]
        if not keyword_set(nolines) then oplot,[break1,break1],[1e-20,1e5],color=orange,line=2
     end
     2: begin 
        pow2=newp[3]
        pow2err=perror[*,3]
        pow3=newp[5]
        pow3err=perror[*,5]
        break1=newp[2]
        break1err=perror[*,2]
        break2=newp[4]
        break2err=perror[*,4]
        
        w1=where(time gt 0 and time lt break1,nw1)
        w2=where(time ge break1 and time lt break2,nw2)
        w3=where(time ge break2)
        t=[time[w1],break1,time[w2],break2,time[w3]]
;        t[0,*]=[tstart[w1],break1,tstart[w2],break2,tstart[w3]]
;        t[1,*]=[tstop[w1],break1,tstop[w2],break2,tstop[w3]]
        y=bkn2pow(t,newp)
        pnames=['Pow1','Breaktime1','Pow2','Breaktime2','Pow3']
        leg=['Pow1'+' = '+sigfig(pow1,3), $
             'Breaktime = '+sigfig(break1,5),$
             'Pow2'+' = '+sigfig(pow2,3),$
             'Breaktime2 = '+sigfig(break2,5),$
             'Pow3'+' = '+sigfig(pow3,3)]
        if not noerr then $
           leg=leg+[' !S!E+'+sigfig(pow1err[1],3)+' !R!I-'+sigfig(pow1err[0],3),$
                    ' !S!E+'+sigfig(break1err[1],5)+' !R!I-'+sigfig(break1err[0],5),$
                    ' !S!E+'+sigfig(pow2err[1],3)+' !R!I-'+sigfig(pow2err[0],3),$
                    ' !S!E+'+sigfig(break2err[1],5)+' !R!I-'+sigfig(break2err[0],5),$
                    ' !S!E+'+sigfig(pow3err[1],3)+' !R!I-'+sigfig(pow3err[0],3) ]
        if not keyword_set(nolines) then begin 
           oplot,[break1,break1],[1e-20,1e5],color=orange,line=2
           oplot,[break2,break2],[1e-20,1e5],color=orange,line=2
        endif 
        
     end
     3: begin 
        pow2=newp[3]
        pow2err=perror[*,3]
        pow3=newp[5]
        pow3err=perror[*,5]
        pow4=newp[7]
        pow4err=perror[*,7]
        break1=newp[2]
        break1err=perror[*,2]
        break2=newp[4]
        break2err=perror[*,4]
        break3=newp[6]
        break3err=perror[*,6]

        w1=where(time ge 0 and time lt break1,nw1)
        w2=where(time ge break1 and time lt break2,nw2)
        w3=where(time ge break2 and time lt break3,nw3)
        w4=where(time ge break3,nw4)
        t=[time[w1],break1,time[w2],break2,time[w3],break3,time[w4]]
;        t[0,*]=[tstart[w1],break1,tstart[w2],break2,tstart[w3],break3,tstart[w4]]
;        t[1,*]=[tstop[w1],break1,tstop[w2],break2,tstop[w3],break3,tstop[w4]]
        y=bkn3pow(t,newp)
        pnames=['Pow1','Breaktime1','Pow2','Breaktime2','Pow3','Breaktime3','Pow4']
        leg=['Pow1'+' = '+sigfig(pow1,3), $
             'Breaktime = '+sigfig(break1,5),$
             'Pow2'+' = '+sigfig(pow2,3),$
             'Breaktime2 = '+sigfig(break2,5),$
             'Pow3'+' = '+sigfig(pow3,3),$
             'Breaktime3 = '+sigfig(break3,5),$
             'Pow4'+' = '+sigfig(pow4,3)]
        if not noerr then $
           leg=leg+[' !S!E+'+sigfig(pow1err[1],3)+' !R!I-'+sigfig(pow1err[0],3),$
                    ' !S!E+'+sigfig(break1err[1],5)+' !R!I-'+sigfig(break1err[0],5),$
                    ' !S!E+'+sigfig(pow2err[1],3)+' !R!I-'+sigfig(pow2err[0],3),$
                    ' !S!E+'+sigfig(break2err[1],5)+' !R!I-'+sigfig(break2err[0],5),$
                    ' !S!E+'+sigfig(pow3err[1],3)+' !R!I-'+sigfig(pow3err[0],3),$
                    ' !S!E+'+sigfig(break3err[1],5)+' !R!I-'+sigfig(break3err[0],5),$
                    ' !S!E+'+sigfig(pow4err[1],3)+' !R!I-'+sigfig(pow4err[0],3)]
        if not keyword_set(nolines) then begin 
           oplot,[break1,break1],[1e-20,1e5],color=orange,line=2
           oplot,[break2,break2],[1e-20,1e5],color=orange,line=2
           oplot,[break3,break3],[1e-20,1e5],color=orange,line=2
        endif 
        
     end
     4: begin 
        pow2=newp[3]
        pow2err=perror[*,3]
        pow3=newp[5]
        pow3err=perror[*,5]
        pow4=newp[7]
        pow4err=perror[*,7]
        pow5=newp[9]
        pow5err=perror[*,9]
        break1=newp[2]
        break1err=perror[*,2]
        break2=newp[4]
        break2err=perror[*,4]
        break3=newp[6]
        break3err=perror[*,6]
        break4=newp[8]
        break4err=perror[*,8]

        w1=where(time ge 0 and time lt break1,nw1)
        w2=where(time ge break1 and time lt break2,nw2)
        w3=where(time ge break2 and time lt break3,nw3)
        w4=where(time ge break3 and time lt break4,nw4)
        w5=where(time ge break4,nw5)
        t=[time[w1],break1,time[w2],break2,time[w3],break3,time[w4],break4,time[w5]]
;        t[0,*]=[tstart[w1],break1,tstart[w2],break2,tstart[w3],break3,tstart[w4]]
;        t[1,*]=[tstop[w1],break1,tstop[w2],break2,tstop[w3],break3,tstop[w4]]
        y=bkn4pow(t,newp)
        pnames=['Pow1','Breaktime1','Pow2','Breaktime2','Pow3','Breaktime3','Pow4','Breaktime4','Pow5']
        leg=['Pow1'+' = '+sigfig(pow1,3), $
             'Breaktime = '+sigfig(break1,5),$
             'Pow2'+' = '+sigfig(pow2,3),$
             'Breaktime2 = '+sigfig(break2,5),$
             'Pow3'+' = '+sigfig(pow3,3),$
             'Breaktime3 = '+sigfig(break3,5),$
             'Pow4'+' = '+sigfig(pow4,3),$
             'Breaktime4 = '+sigfig(break4,5),$
             'Pow5'+' = '+sigfig(pow5,3)]
        if not noerr then $
           leg=leg+[' !S!E+'+sigfig(pow1err[1],3)+' !R!I-'+sigfig(pow1err[0],3),$
                    ' !S!E+'+sigfig(break1err[1],5)+' !R!I-'+sigfig(break1err[0],5),$
                    ' !S!E+'+sigfig(pow2err[1],3)+' !R!I-'+sigfig(pow2err[0],3),$
                    ' !S!E+'+sigfig(break2err[1],5)+' !R!I-'+sigfig(break2err[0],5),$
                    ' !S!E+'+sigfig(pow3err[1],3)+' !R!I-'+sigfig(pow3err[0],3),$
                    ' !S!E+'+sigfig(break3err[1],5)+' !R!I-'+sigfig(break3err[0],5),$
                    ' !S!E+'+sigfig(pow4err[1],3)+' !R!I-'+sigfig(pow4err[0],3),$
                    ' !S!E+'+sigfig(break4err[1],5)+' !R!I-'+sigfig(break4err[0],5),$
                    ' !S!E+'+sigfig(pow5err[1],3)+' !R!I-'+sigfig(pow5err[0],3)]
        if not keyword_set(nolines) then begin 
           oplot,[break1,break1],[1e-20,1e5],color=orange,line=2
           oplot,[break2,break2],[1e-20,1e5],color=orange,line=2
           oplot,[break3,break3],[1e-20,1e5],color=orange,line=2
           oplot,[break4,break4],[1e-20,1e5],color=orange,line=2
        endif 
        
     end
  endcase 
  
  if norm gt 100 then sci=1 else sci=0
  normleg='Norm = '+sigfig(norm,3,sci=sci)
  if not noerr then $
     normleg=normleg+' !S!E+'+sigfig(normerr[1],3,sci=sci)+' !R!I-'+sigfig(normerr[0],3,sci=sci)
  
  leg=[leg,$
       normleg,$
       sym.chi+'!U2!N/dof = '+sigfig(chisq/dof,4),$
       'dof = '+ntostr(fix(dof))]
  
;  if not keyword_set(nolines) then 
  oplot,t,y,color=green,thick=1
  if not keyword_set(noleg) then legend,leg,box=0,/top,/right,charsize=charsize


  return
end 

pro fit_lc_sub,file,newp=newp,yfit=yfit,t=t,perror=perror,lc=lc,oldfile=oldfile,name=name,phil=phil,_extra=_extra,arrowsize=arrowsize,xrange=xrange,nohard=nohard,xtitle=xtitle,qdp=qdp,int=int,flux=flux,nocolor=nocolor
  
  if n_elements(lc) eq 0 then lc=lcout2fits(file,qdp=qdp)
  time=lc.time
  tstarted=lc.tstart
  tstoped=lc.tstop
  cts=lc.src_rate*flux
  err=lc.src_rate_err*flux
  type=lc.type
  bg=lc.tot_back_cts
  src=lc.src_counts
  sig=lc.det_sig
  expt=lc.exptime
;  sigma=src/sqrt(src+bg*2)
  sigma=sig
  chisq2=0.
  dof2=0.
  
  ul=where(err le 0,nul)
  wdet=where(err gt 0)
  
  timerr=fltarr(2,n_elements(time))
  timerr[0,*]=time-tstarted
  timerr[1,*]=tstoped-time
;  timerr=((time-tstarted)+(tstoped-time))/2.
  type=fix(type)
  w=where(cts gt 0 and finite(err) and err ne 0)
  time=time[w] & timerr=timerr[*,w] & cts=cts[w] & err=err[w] & type=type[w]
  back=bg                       ;*expt

;  if cts[1]*expt[1] lt 12 then cash=1 else cash=0
;  if src[1] lt 12 then cash=1 else cash=0
;  if cash then begin
;     print,'USING CSTATS IN XSPEC'
;     fit_lc_xspec,newp,perror,chisq,dof,yfit,slope=1,/cash,/noplot
;     status=1
;     np=n_elements(newp)
;     case np of
;        2: breaks=0
;        4: breaks=1
;        6: breaks=2
;        8: breaks=3
;     endcase 
;  endif else begin

  refit:
;     plot_base_lc,_extra=_extra
  
;  multiplot2,/reset,/default
  erase
  if nohard eq 0 then begin 
     multiplot2,[1,2],/init
     multiplot2
  endif 
  plot_like_qdp,_extra=_extra,lc=lc[wdet],title=name,arrowsize=arrowsize,xrange=xrange,xtitle=xtitle,noxaxis=noxaxis,qdp=qdp,file=file,flux=flux,yrange=yrange,nocolor=nocolor

  bt=0d
  !mouse.button=0

  while (!MOUSE.button NE 4) do begin
     
     print,'Click on estimate breaktime, or right click to continue'
     cursor,xxx,yyy,/wait,/change
     if !mouse.button ne 4 then begin 
        oplot,[xxx,xxx],[1e-20,1e4],color=!orange
        print,round(xxx)
        bt=[bt,xxx]
     endif
  endwhile
  breaks=n_elements(bt)-1
  if breaks gt 0 then bt=round(bt[1:*])

  slope=0d
  mintime=ntostr(min(time))
  maxtime=ntostr(max(time))
  for i=0,breaks do begin 
     if i lt breaks then begin
        if i ne 0 then time1=ntostr(bt[i-1]) else time1=mintime
        wtime2=ntostr(bt[i])
        if i lt breaks-1 then time2=ntostr(bt[i+1]) else time2=maxtime
     endif else begin 
        if breaks gt 0 then time1=ntostr(bt[i-1]) else time1=mintime
        wtime2=maxtime
     endelse 
     w=where(time ge time1-10. and time le wtime2+10. and err ne 0,nw)
     mw=max(w)+1
     w2=w
     if mw lt n_elements(time) then $
        w2=[w,mw]
     f=linfit(alog10(time[w]),alog10(cts[w]))
     oplot,time[w2],10^f[0]*time[w2]^f[1]
     print,time1,' ',wtime2
     sl=-f[1]
     slope=[slope,sl]
;      minslope=ntostr(sl-2)
;      maxslope=ntostr(sl+2)
     if i eq 0 then norm=10^f[0]*1^f[1]
  endfor 
  slope=slope[1:*]

  if not finite(norm) then norm=max(cts)

;  norm=1e31;max(cts)*10.
  case breaks of
     0: begin
;           mo='pow'
        mo='intpow'
        p=[norm,slope]
        pnames=['norm','pow']
        pmin0=[0.,-10.]
        log=[1,0]
     end 
     1: begin 
;           mo='bknpow'
        mo='intbknpow'
        p=[norm,slope[0],bt,slope[1]]
        pnames=['norm','pow1','break','pow2']
        pmin0=[0.,-10.,0,-10.]
        log=[1,0,1,0]
     end 
     2: begin 
;           mo='bkn2pow'
        mo='intbkn2pow'
        p=[norm,slope[0],bt[0],slope[1],bt[1],slope[2]]
        pnames=['norm','pow1','break1','pow2','break2','pow3']
        pmin0=[0.,-10.,0,-10.,0.,-10.]
        log=[1,0,1,0,1,0]
     end 
     3: begin 
;           mo='bkn3pow'
        mo='intbkn3pow'
        p=[norm,slope[0],bt[0],slope[1],bt[1],slope[2],bt[2],slope[3]]
        pnames=['norm','pow1','break1','pow2','break2','pow3','break3','pow4']
        pmin0=[0.,-10.,0,-10.,0.,-10.,0.,-10.]
        log=[1,0,1,0,1,0,1,0]
     end 
     4: begin 
;           mo='bkn4pow'
        mo='intbkn4pow'
        p=[norm,slope[0],bt[0],slope[1],bt[1],slope[2],bt[2],slope[3],bt[3],slope[4]]
        pnames=['norm','pow1','break1','pow2','break2','pow3','break3','pow4','break4','pow5']
        pmin0=[0.,-10.,0,-10.,0.,-10.,0.,-10.,0.,-10.]
        log=[1,0,1,0,1,0,1,0,1,0]
     end 
  endcase
  
  fit_pow_model,time[wdet],cts[wdet],timerr[*,wdet],err[wdet],p,mo,pnames,yfit,newp,perror,chisq,dof,weights,src[wdet],back[wdet],status=status
;  endelse 
  
  if status eq 0 then goto,refit
  
  perror0=dblarr(2,n_elements(perror))
  perror0[0,*]=perror
  perror0[1,*]=perror

  plot_lcfit_results,lc,newp,perror0,chisq,dof,breaks,leg,pnames,/noerr,noleg=noleg,charsize=charsize,ps=ps,nolines=nolines,name=name,flux=flux,nocolor=nocolor
  if nohard then nwhard=0 else whard=where(oldlc.tot_hard lt 1e3,nwhard)
;  whard=where(lc.tot_hard lt 1e3,nwhard)
  if nwhard gt 0 and not keyword_set(nohard) then begin 
     multiplot2
     hyrange=[max([min(lc[whard].tot_hard-lc[whard].tot_hard_err),1e-2]),max(lc[whard].tot_hard+lc[whard].tot_hard_err)]
     ploterror,lc[whard].time,lc[whard].tot_hard,lc[whard].tot_hard_err,psym=3,/nohat,xtitle='Time since BAT trigger (s)',ytitle='hardness ratio',charsize=1.,/xlog,/ylog,yrange=hyrange,xrange=xrange
     for r=0,n_elements(whard)-1 do oplot,[lc[whard[r]].tstart,lc[whard[r]].tstop],[lc[whard[r]].tot_hard,lc[whard[r]].tot_hard]  
  endif 
;  multiplot2,/reset

  ;;;F-test
  print,dof
  if chisq2 ne 0 then begin
     num=n_elements(lc)
     m=num-dof
     ddof=dof2-dof
     f=ftest(chisq2,chisq,m,num,ddof)
     print,'F-test:  ',f
;     stop
  endif 
  chisq2=chisq
  dof2=dof
  
  doagain='y'
  input,'Is this an acceptable fit? (y/n) ',doagain,'y'
  if doagain eq 's' then stop
  if doagain eq 'n' then begin
     slope=-999
     goto,refit
  endif 
  
  perror0=perror
  print
  print,'Calculating 90% confidence errors'
  tt=dblarr(2,n_elements(time))
  tt[0,*]=time-timerr[0,*]
  tt[1,*]=time+timerr[1,*]
  
  multiplot2
  multiplot2,/reset
  erase
  !p.multi=0

  if finite(perror[0]) then begin 
     multiplot2,/default
;     delchi0=chisqr_cvf(0.1,n_elements(newp)-1)
;     if n_elements(newp) eq 2 then delchi0=chisqr_cvf(0.1,1)
;     if n_elements(newp) ge 4 then delchi0=chisqr_cvf(0.1,3)
;     print,delchi0
;     conf_error,tt,cts,err,newp,perror0,mo,perror2,bestfit,pvarunit,bestchisq,yfit,log=log,pmin0=pmin0,delchi0=delchi0,/doplot
     print,'Using Monte Carlo error method with 1000 simulations'
     
     lc_monte_pow,lc,newp,['norm',pnames],chisq,dof,perror,/noplot,ps=ps,/nowrite,nsim=1000,int=int
     
     k=get_kbrd(10)
;  colprint,newp,perror0,perror2[0,*],perror2[1,*],bestfit,pvarunit,bestchisq
;  newp=bestfit*1.
;  for i=0,n_elements(newp)-1 do perror[i]=mean(perror2[*,i])
;     perror=perror2
  endif else perror=dblarr(2,n_elements(perror))
;  chisq=total((cts-yfit)^2./err^2)
;  oplot,time,yfit,color=!red
  newp=newp*1d
  
;  erase  
;  plot_base_lc
  
  erase
  if not nohard then begin 
     multiplot2,[1,2],/init
     multiplot2
  endif 
  plot_like_qdp,file=oldfile,title=name,_extra=_extra,xrange=xrange,xtitle=xtitle,noxaxis=noxaxis,qdp=qdp,flux=flux,yrange=yrange,nocolor=nocolor
  plot_lcfit_results,lc,newp,perror,chisq,dof,breaks,leg,pnames,noleg=noleg,charsize=charsize,ps=ps,nolines=nolines,_extra=_extra,name=name,flux=flux,nocolor=nocolor
  if nohard then nwhard=0 else whard=where(oldlc.tot_hard lt 1e3,nwhard)
;  whard=where(lc.tot_hard lt 1e3,nwhard)
  if nwhard gt 0 then begin 
     multiplot2                 ;,ydowngap=0.08
     hyrange=[max([min(lc[whard].tot_hard-lc[whard].tot_hard_err),1e-2]),max(lc[whard].tot_hard+lc[whard].tot_hard_err)]
     ploterror,lc[whard].time,lc[whard].tot_hard,lc[whard].tot_hard_err,psym=3,/nohat,xtitle='Time since BAT trigger (s)',ytitle='hardness ratio',charsize=1.,/xlog,/ylog,yrange=hyrange,xrange=xrange
     for r=0,n_elements(whard)-1 do oplot,[lc[whard[r]].tstart,lc[whard[r]].tstop],[lc[whard[r]].tot_hard,lc[whard[r]].tot_hard]
  endif 
  multiplot2,/reset
  !p.multi=0
  t=time
  ;;write output ps file
  begplot,name='lc_fit_plot.ps',/landscape,/color
  plot_like_qdp,file=oldfile,lc=oldlc,title=name,_extra=_extra,pmulti=pmulti,symsize=symsize,xrange=xrange,xtitle=xtitle,noxaxis=noxaxis,qdp=qdp,flux=flux,yrange=yrange,nocolor=nocolor
;  plot_like_qdp,_extra=_extra,name=name,phil=phil
;  plot_lcfit_results,lc,newp,perror,chisq,dof,breaks,leg,pnames,noleg=noleg,charsize=charsize
  plot_lcfit_results,lc,newp,perror,chisq,dof,breaks,leg,pnames,_extra=_extra,noleg=noleg,charsize=charsize,ps=ps,nolines=nolines,name=name,flux=flux,nocolor=nocolor
;  oplot,time,yfit,color=!green
;  legend,leg,box=0,/top,/right
  endplot
  
  ;;write output fit file
  print,'writing out lc_fit_out_idl_int'+int+'.dat'
  openw,lun,'lc_fit_out_idl_int'+int+'.dat',/get_lun
  norm=newp[0]
  normerr=perror[*,0]
  for i=0,n_elements(newp)-2 do begin
     j=i+1
     printf,lun,pnames[i]+' '+ntostr(newp[j])+' '+ntostr(perror[0,j])+' '+ntostr(perror[1,j])
  endfor
  printf,lun,'Norm '+ntostr(norm)+' '+ntostr(normerr[0])+' '+ntostr(normerr[1])
  printf,lun,'Chisq '+ntostr(chisq)
  printf,lun,'dof '+ntostr(dof)
  close,lun
  free_lun,lun
  
  return
end 
pro fit_lc_wrapper,phil=phil,qdp=qdp
  
  cd,!mdata
  g=0
  if n_elements(dir) eq 0 then dir=file_search('GRB*')
  nw=n_elements(dir)

  for i=g,nw-1 do begin 
     cd,dir[i]
     print
     print,dir[i],i
     if exist('lc_newout.txt') or keyword_set(phil) then $
        fit_lc,phil=phil,qdp=qdp,name=dir[i]+' '+ntostr(i)
     cd,'..'
  endfor 
  
  return
end 
pro fit_lc,file=file,name=name,newp=newp,perror=perror,_extra=_extra,phil=phil,pmulti=pmulti,lc=lc,noleg=noleg,ps=ps,noxaxis=noxaxis,ytitle=ytitle,nocolor=nocolor,justplot=justplot,noinit=noinit,nohard=nohard,xtitle=xtitle,qdp=qdp,nolines=nolines,int=int,flux=flux
  
  if not keyword_set(flux) then flux=1.
  red=!p.color
  if not keyword_set(nocolor) then red=!red
  usefile=0
  if n_elements(file) eq 0 then phil=1 else usefile=1
  simpctable
  defsymbols
;  !p.multi=0
  if n_elements(file) eq 0 then file='lc_newout.txt'
  if keyword_set(phil) then begin
     convert_phil2lcout,qdp=qdp
     file='lc_newout_phil.txt'
     if exist('lc_newout_chandra.txt') then begin
        spawn,'cat lc_newout_phil.txt lc_newout_chandra.txt > lc_newout_phil2.txt'
        spawn,'cp lc_newout_phil2.txt lc_newout_phil.txt'
     endif 
     nohard=1
;     phil=0
  endif 
  oldlc=lcout2fits(file,qdp=qdp)
;  lc=lcout2fits(file)
  lc=oldlc
  xrange=[lc[0].tstart,lc[n_elements(lc)-1].tstop]

  wdet=where(lc.src_rate_err gt 0)
  if not keyword_set(nohard) then begin
     xtitle=' '
     nohard=0
  endif else xtitle='Time since BAT trigger (s)' 
     
  ep=0
  if keyword_set(ps) then begin
     begplot,name='lc_fit_plot.ps',/color,/land
     ep=1
  endif 
  
  if exist('lc_newout.txt') or keyword_set(phil) or keyword_set(qdp) or usefile then begin 
     slope=-999
     nfl=''
     if n_elements(pmulti) eq 0 then begin 
        if not keyword_set(noinit) then begin 
           erase
           if not nohard then  multiplot2,[1,2],/init
           noinit=0
           ydowngap=0
        endif else ydowngap=0.08
        if not nohard then multiplot2
        multi=1
     endif else begin
        multi=0
        noinit=0
     endelse 

     plot_like_qdp,file=file,title=name,lc=oldlc,_extra=_extra,pmulti=pmulti,symsize=symsize,ytitle=ytitle,nocolor=nocolor,xrange=xrange,noxaxis=noxaxis,xtitle=xtitle,xtickname=xtickname,qdp=qdp,flux=flux,yrange=yrange
     ul=where(lc.src_rate_err le 0,nul)
     if nul gt 0 then begin
        plotsym,1,3,thick=3
        plots,lc[ul].time,lc[ul].src_rate*flux,psym=8,color=red
        for uu=0,nul-1 do oplot,[lc[ul[uu]].tstart,lc[ul[uu]].tstop],[lc[ul[uu]].src_rate,lc[ul[uu]].src_rate]*flux,color=red
     endif 
        
;     replot_xrt_lc,time,timerr,cts,err,file='lc_newout.txt',title=dir[i]
     if n_elements(int) eq 0 then begin 
        int='7'
        ffile='lc_fit_out_idl_int7.dat'
        if not exist(ffile) then ffile='lc_fit_out_idl_int6.dat'
        if not exist(ffile) then ffile='lc_fit_out_idl_int5.dat'
        if not exist(ffile) then ffile='lc_fit_out_idl_int4.dat'
        if not exist(ffile) then ffile='lc_fit_out_idl_int3.dat'
     endif else ffile='lc_fit_out_idl_int'+int+'.dat'
     go=0
     if exist('lc_newout_noflares.txt') and not exist('flares_gtis.dat') then go=1
     if exist(ffile) and not go then begin 
        read_lcfit,ffile,pname,newp,perror,chisq,dof,breaks
        if pname[0] ne 'nofit' then begin 
           plot_lcfit_results,lc,newp,perror,chisq,dof,breaks,leg,pnames,_extra=_extra,noleg=noleg,charsize=charsize,ps=ps,nocolor=nocolor,nolines=nolines,name=name,flux=flux
           if multi and not nohard  then multiplot2,ydowngap=ydowngap
           if nohard then nwhard=0 else whard=where(oldlc.tot_hard lt 1e3,nwhard)
           if nwhard gt 0 then begin 
              hyrange=[max([min(oldlc[whard].tot_hard-oldlc[whard].tot_hard_err),1e-2]),max(oldlc[whard].tot_hard+oldlc[whard].tot_hard_err)]
              ploterror,oldlc[whard].time,oldlc[whard].tot_hard,oldlc[whard].tot_hard_err,psym=3,/nohat,xtitle='Time since BAT trigger (s)',ytitle='hardness ratio',charsize=1.,/xlog,/ylog,yrange=hyrange,xrange=xrange,xtickf='loglabels';,xtickname=xtickname
              for r=0,n_elements(whard)-1 do oplot,[oldlc[whard[r]].tstart,oldlc[whard[r]].tstop],[oldlc[whard[r]].tot_hard,oldlc[whard[r]].tot_hard]
           endif 
           if multi and noinit eq 0 and not nohard then multiplot2,/reset
           if keyword_set(justplot) then begin 
              if keyword_set(ps) then endplot
              return
           endif 
           print,'Use previous fit? (y/n)'
           pfit='y'
           pfit=get_kbrd(10)
           if pfit eq 's' then stop
        endif else begin 
           pfit='n'
           nfl='p'
        endelse 
     endif else pfit='n'
     if pfit eq 'n' then begin 
        
        newfile='lc_newout_noflares.txt'
        ans='n'
        if exist(newfile) then begin
           print,'Use existing noflares filter? (y/n)'
           ans=get_kbrd(10)
        endif
        if ans eq 'y' then lc=lcout2fits(newfile,qdp=qdp)
        if ans eq 'n' then begin 
           print,'remove flares or p to skip? (y/n/p)'
           fl=get_kbrd(10)
           fit_noflares_again:
           if fl eq 'y' then begin
              if multi then begin 
                 erase
                 if not nohard then begin 
                    multiplot2,[1,2],/init
                    multiplot2
                    multiplot2
                 endif 
              endif 
              if nohard then nwhard=0 else whard=where(oldlc.tot_hard lt 1e3,nwhard)
              if nwhard gt 0 then begin 
                 hyrange=[max([min(oldlc[whard].tot_hard-oldlc[whard].tot_hard_err),1e-2]),max(oldlc[whard].tot_hard+oldlc[whard].tot_hard_err)]
                 ploterror,oldlc[whard].time,oldlc[whard].tot_hard,oldlc[whard].tot_hard_err,psym=3,/nohat,xtitle='Time since BAT trigger (s)',ytitle='hardness ratio',charsize=1.,/xlog,/ylog,yrange=hyrange,xrange=xrange
                 for r=0,n_elements(whard)-1 do oplot,[oldlc[whard[r]].tstart,oldlc[whard[r]].tstop],[oldlc[whard[r]].tot_hard,oldlc[whard[r]].tot_hard]
              endif 
              if multi and not nohard then begin 
                 multiplot2,/reset
                 multiplot2,[1,2],/init
                 multiplot2
              endif 
              fit_noflares2,file,/small,lc=lc,qdp=qdp,_extra=_extra,title=name
           endif else begin 
              if exist(newfile) and newfile eq 'lc_newout_noflares.txt' then begin 
                 spawn,'rm '+newfile
                 spawn,'rm flares_gtis.dat'
              endif 
;              newfile='lc_newout.txt'
              newfile=file
;              lc=oldlc
           endelse 
           if fl eq 's' then stop
        endif else fl=''
        if ans ne 'p' and fl ne 'p' then begin 
           erase
           if multi and not nohard then begin 
              multiplot2,[1,2],/init
              multiplot2
           endif 
           plot_like_qdp,file=newfile,lc=lc,_extra=_extra,title=name,symsize=symsize,xrange=xrange,xtitle=xtitle,noxaxis=noxaxis,qdp=qdp,flux=flux,yrange=yrange,nocolor=nocolor

           if multi and not nohard then multiplot2
           if nohard then nwhard=0 else whard=where(oldlc.tot_hard lt 1e3,nwhard)
           if nwhard gt 0 then begin 
              hyrange=[max([min(oldlc[whard].tot_hard-oldlc[whard].tot_hard_err),1e-2]),max(oldlc[whard].tot_hard+oldlc[whard].tot_hard_err)]
              ploterror,oldlc[whard].time,oldlc[whard].tot_hard,oldlc[whard].tot_hard_err,psym=3,/nohat,xtitle='Time since BAT trigger (s)',ytitle='hardness ratio',charsize=1.,/xlog,/ylog,yrange=hyrange,xrange=xrange
              for r=0,n_elements(whard)-1 do oplot,[oldlc[whard[r]].tstart,oldlc[whard[r]].tstop],[oldlc[whard[r]].tot_hard,oldlc[whard[r]].tot_hard]
           endif 
           if multi and not nohard then multiplot2,/reset
           print,'Flare removal (or not) is sufficient? (y/n)'
           ag=get_kbrd(10)
           if ag eq 'n' then begin
              fl='y'
              goto,fit_noflares_again
           endif 
           if multi and not nohard then multiplot2,/reset
           if ag eq 's' then stop
        endif
        if fl[0] ne 'n' and fl[0] ne 'y' and nfl eq 'p' then fl='p'
;        print,'type to continue (s to stop)'
;        k=get_kbrd(10)
;        if k eq 's' then stop
        if fl ne 'p' then begin
           print,newfile
           fit_lc_sub,newfile,newp=newp,lc=lc,xrange=xrange,nohard=nohard,qdp=qdp,_extra=_extra,name=name,oldfile=file,int=int,flux=flux,nocolor=nocolor
           print,'type to continue'
           k=get_kbrd(10)
           if k eq 's' then stop
        endif else begin
           print,'writing out lc_fit_out_idl_int'+int+'.dat'
           openw,lun,'lc_fit_out_idl_int'+int+'.dat',/get_lun
           printf,lun,'no fit'
           close,lun
           free_lun,lun
        endelse 
     endif else begin
        print,'type to continue'
        k=get_kbrd(10)
     endelse
;        plot_like_qdp,t,yfit,dir[i]
  endif  
  
  if multi and noinit eq 0 then multiplot2,/default
  return  
end

; defaults write com.apple.x11 wm_click_through -bool true 
