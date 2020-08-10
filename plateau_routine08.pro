pro plateau_routine08
;  cleanflares                 
  fit_lc,/phil  ;LC fit
  read_lcfit,'lc_fit_out_idl_int7.dat',pname,p,perror 
  format_error,pname,p,perror
  print,'type to continue'
  k=get_kbrd(10)
  norm_plateau_end=p[0]*p[2]^(p[1]) 
  plateaufit,norm_plateau_end,p[3],p[2]
  print,'type to continue'
  k=get_kbrd(10)
  fit_lc,/phil
  readcol,'lc_fit_out_idl_int7.dat',pname,p,format='(a,d)'
  print,p[1],p[0],p[2],p[3],format='(d10.1,d10.2,d7.2,d6.0)'
  spawn,'mv PCCURVE.qdp.temp PCCURVE.qdp'
  spawn,'mv WTCURVE.qdp.temp WTCURVE.qdp'
  spawn,'mv lc_fit_out_idl_int7.dat_temp lc_fit_out_idl_int7.dat'
end

