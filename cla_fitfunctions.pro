function fit_offpow,t,p,perror
  ;Simple version, where I don't specify the errors
  norm=p[0]
  pow1=p[1]
  T0=p[2]
  yfit=norm*(t+T0)^(-pow1)
  if n_elements(perror) gt 0 then begin
    normerr=perror[0]
    powerr=perror[1]
    T0err=perror[2]
    dydn=(t+T0)^(-pow1)
    dydTO=-norm*pow1*(t+T0)^(-pow1-1)
    yerr2=normerr^2*dydn^2+T0err^2*dyT0^2
    yerror=sqrt(yerr2)
 endif
  return,yfit
end



