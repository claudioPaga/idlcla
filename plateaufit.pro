@cla_fitfunctions
pro plateaufit,norm,pow1,tplateau
  if n_params() lt 3 then begin
     print,'Please enter correct parameters: plateaufit,norm,pow_index,t_plateau'
     return
  endif
  ;Here I get the GRB name'
  spawn,'pwd > nomedir.txt'
  readcol,'nomedir.txt',nomedir,format='(a)'
  grbname=strsplit(nomedir,'/',/extract)     
  grbname=strtrim(grbname[4],2)
    
                                ;pow1=Decay index of post-plateau phase from original LC fit
    
  readcol,'lc_priorfit_noflares.txt',t,tstart,tstop,cts,err,format='(d,d,d,d,d)'
                                ;sometimes the last point in the LC is an upper limit
                                ;and its error is set to 0 by Phil
                                ;mpfitfun does not like errors=0, so bnelow I remove those lines
  
  index=where(err ne 0.)
  t=t[index]
  tstart=tstart[index]
  tstop=tstop[index]
  cts=cts[index]
  err=err[index]
  
  terr=tstop-tstart
  p=[norm,pow1,tplateau];this is the parameters initial guess
  
  ;Here I format the parameters.  In particular, pow1, the decay index after the plateau phase is actually a fixed par
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                       limits:[0.D,0]}, 3)
  
  parinfo[0].limited[0] = 1
  parinfo[0].limits[0]  = 0.D ;Here I say that the normalization must be greater than 0
  parinfo[1].fixed = 1
  parinfo[2].limited[0] = 1
  parinfo[2].limits[0]  = 0.D ;Here I say that the prior time should be greater than 1000.
  parinfo[*].value = p
  
;  A total of 3 parameters, with starting values of pow1,
;  1, and 1000 are given.  The first parameter
;  is fixed at a value of pow1, and the last parameter is
;  constrained to be above 1000.
;  stop 
  ecco=mpfitfun('fit_offpow',t,cts,err,parinfo=parinfo,perror=perror,bestnorm=chisq,dof=dof)  
 ; ecco=mpfitfun('fit_offpow',t,cts,err,p) 
  print,'Norm,PI,T0: ',ecco
  print,'Param_errors: ',perror
  print,'Corrected Param errors: ',perror*sqrt(chisq/dof)
  ;Here I look for the file with times within 45 degrees of FoV, to get the predicted countrates at these times
  filesearch='/home/pagani/xrt/plateau/asflown/'+grbname+'_times_fov.txt'
  readcol,filesearch,ra,dec,distance,fovtimesmet,fovtimestart,fovtimeend,skipline=2,format='(d,d,d,d,d,d)'
  withintimes=where(fovtimeend lt ecco[2])
  if withintimes[0] ne -1 then begin 
     fovtimes_start=fovtimestart(withintimes)
     fovtimes_end=fovtimeend(withintimes)
     fovpredict=dblarr(n_elements(fovtimes_start),2)
     
     ;Here I check if the prior time is inside one of the snapshots (that is, if the GRB was inside the FOV at the prior time)
     
     checcafov=0
     for k=0,n_elements(fovtimes_start)-1 do begin
        if ecco[2] lt fovtimes_start(k) and ecco[2] gt fovtimes_end(k) then checcafov=1
     endfor
     if checcafov eq 1 then print,'YES!!!!  GRB inside the BAT FoV at the prior time' else print,'NO!!!!! GRB NOT inside BAT FoV at prior time'
     
     
     
     for i=0,n_elements(fovtimes_start)-1 do begin
        if ecco[2]-fovtimes_start(i) lt 0 then begin
           fovpredict(i,0)=ecco[0];Predicted CR at start of IN FoV snapshot
           fovpredict(i,1)=ecco[0]*(ecco[2]-fovtimes_end(i))^(-ecco[1]) ;Predicted CR at end of IN FoV snapshot
           if i eq 0 then  print,'CR at T_prior+',ecco[2]-ecco[2],fovpredict(i,0) else print,'CR at T_prior+',ecco[2],fovpredict(i,0)         
        endif else begin
           fovpredict(i,0)=ecco[0]*(ecco[2]-fovtimes_start(i))^(-ecco[1]) ;Predicted CR at start of IN FoV snapshot
           fovpredict(i,1)=ecco[0]*(ecco[2]-fovtimes_end(i))^(-ecco[1]) ;Predicted CR at end of IN FoV snapshot
           print,'CR at T_prior+',ecco[2]-fovtimes_start(i),fovpredict(i,0)
        endelse   
     endfor
     if checcafov eq 1 then fovtimes_start[0]=ecco[2]
     filewritename=grbname+'_cr_extrap.txt'
     openw,filewrite,filewritename,/get_lun
     cr_table=[transpose(fovtimes_start),transpose(fovtimes_end),transpose(fovpredict)]
     printf,filewrite,'FoV_start,  FoV_end,  CR_start,  CR_end'
     printf,filewrite,cr_table
     free_lun,filewrite
     
     
     
  endif else print,'NO TIMES OF GRB IN BAT FOV AFTER T_PRIOR'
  print,'CR_prior  Err_CR_prior T_prior Error_T_prior'
  sp='  '
  print,strtrim(string(ecco[0]),2),sp,strtrim(string(perror[0]*sqrt(chisq/dof)),2),sp,strtrim(string(ecco[2]),2),sp,strtrim(string(perror[2]*sqrt(chisq/dof)),2)
  refit_prior,ecco[2]
  
end
