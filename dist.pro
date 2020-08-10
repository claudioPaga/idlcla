pro dist,redshift,fluence,pi_gamma,break_time
  
  if n_params() lt 3 then begin
     print,''
     print,'WARNING!!!!!   Wrong number of parameters!!!!    To call the procedure use:'
     print,''
     print,'dist,redshift,fluence,bat_PhotonIndex,Jet_break_time'
     print,''
     print,'Please enter Fluence in 10^-6ergs/cm2 units,pi_gamma is the BAT photon index, so its the value from the fit not the PI-1 value and break_time is in seconds, the script converts it into days for jet_break angle calculations'
     return
  endif
  
PI=3.14159
print,''
dL=lumdist(redshift,h0=71,omega_m=0.27,lambda0=0.73)

;Calculate Lum_Distance with units of 10^24 cm

dLcm=dL*3.085677581
print,'dL=',dL,' Mpc','    dL=',dLcm/10000.0,'x10^28cm'
print,''

;Calculate isotropic energy

En_15_150=4*PI*fluence*dLcm*dLcm
print,'E(iso,BAT) no K-corrected: ',En_15_150/1.0e+10,'x10^52 erg'

;Apply simple K-correction, assuming a simple power law spectrum for
;the GRB prompt emission

En_15_150_kcorr=En_15_150*(1+redshift)^(pi_gamma-1)/1.0e+10
print,'E(iso,BAT), K-corrected= ',En_15_150_kcorr,'x10^52 erg'

;Convert break time in days units.  Define param for Jet angle formula

b_time_days=break_time/86400.0
efficency=0.2
density=10.0

angle_jet=0.161*(b_time_days/(1+redshift))^(3.0/8.0)*(efficency*density/En_15_150_kcorr)^(1.0/8.0)
;print,angle_jet
print,''
ang=strtrim(string(angle_jet),2)
deg=strtrim(string(angle_jet/(2*pi)*360.0),2)
print,'Jet-break angle= ',ang,'   (',deg,' degrees)'

;Calculate "true" GRB energy
print,''
E_gamma=(1-cos(angle_jet))*En_15_150_kcorr*100.0
print,'True Energy(BAT)= ',E_gamma,'x10^50 erg'

end
