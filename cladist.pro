function cladist,ra1,dec1,ra2,dec2

  if n_params() lt 4 then begin
     print,'Please input: cladist(ra1,dec1,ra2,dec2) in decimal degrees'
     return,0
 endif

 ra1=ra1*!PI/180.0
 ra2=ra2*!PI/180.0
 dec1=dec1*!PI/180.0
 dec2=dec2*!PI/180.0
 
 d=acos(cos(!PI/2.0-dec1)*cos(!PI/2.0-dec2)+sin(!PI/2.0-dec1)*sin(!PI/2.0-dec2)*cos(ra1-ra2))
 
;calculated d in arcsec 

 jd=d*180.0/!PI*3600.0

  print,'Separation in arcsec: ',jd

  return,jd
end
