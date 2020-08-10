pro refit_prior,tprior
  qdpfiles=findfile('*CURVE.qdp')
  for k=0,n_elements(qdpfiles)-1 do begin
     comcopia='cp '+qdpfiles(k)+' '+qdpfiles(k)+'.temp'
     spawn,comcopia
     readcol,qdpfiles(k),time,tposerr,tnegerr,rate,err,fracexp,bgrate,bgerr,corr_fact,Cts_in_src,bg_in_src,exposure,sigma,format='(d,d,d,d,d,d,d,d,d,d,d,d,d)'
     time=time+tprior
     writecol,qdpfiles(k),time,tposerr,tnegerr,rate,err,fracexp,bgrate,bgerr,corr_fact,Cts_in_src,bg_in_src,exposure,sigma
  endfor
  spawn,'mv lc_fit_out_idl_int7.dat lc_fit_out_idl_int7.dat_temp'
end
