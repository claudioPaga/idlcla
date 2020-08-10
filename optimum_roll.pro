pro optimum_roll

rangc=186.454  
decngc=33.547
openw,allow,'allowed_roll.txt',/get_lun
for i=0,365 do begin
   jd=2454566.5+i               ;This is the Julian date as of April 10th, 2008
   sunpos, jd, rasun, decsun
 ;   # define the Sun and target vector    
   vector,decsun,rasun,vSun
   vector,decngc,rangc,vT
;print,'vSun:',vSun
;print,'vT:',vT
  ;  # get cross product of Sun vector and target vector (vSun x vT) :  result is vector Y

    vY = crossp(vT, vSun)
    vnY=vY/norm(vY)
;print,'vnY',vnY
   ; # get cross product of normalized Y and target vector: result is vector Z

    vZ=  crossp(vnY, vT)
    vnZ=vZ/norm(vZ)
;print,'vnZ',vnZ
 ;   # newroll = -atan(vnY[2]/vnZ[2]), to get proper newroll use function atan2
 ;   # JAK Note: For some reason this doesnt work - so I got rid of the sign and it did.

    if (vnY[2] ne 0) and (vnZ[2] ne 0) then begin
        newroll = atan(vnY[2],vnZ[2])
        newroll = newroll/!dtor
     endif else begin
;        # roll is not uniquely defined, arbitrarily pick 0.0 or 180.0
        newroll = 0
        if (vSun[0]*(-cos(newroll)*sin(decngc*!dtor)*cos(rangc*!dtor) - sin(newroll)*sin(rangc*!dtor)) + vSun[1]*(-cos(newroll)*sin(decngc*!dtor)*sin(rangc*!dtor) + sin(newroll)*cos(rangc*!dtor)) + vSun[2]*cos(newroll)*cos(decngc*!dtor) < 0.0) then newroll=180.0
  ;          # 0.0 would put the sun vector in -Z, so use 180.0 deg
     endelse 
    if (newroll lt 0.0) then newroll = newroll+360.0
    ;return newroll*!dtor
if i lt 265 then begin
printf,allow,'Doy',i+100,'  Optimal roll:',newroll
endif else begin
   printf,allow,'Doy',i-265,'  Optimal roll:',newroll
                                
endelse
endfor
free_lun,allow
end
