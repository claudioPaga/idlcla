pro roll_check,year,doy,ra_ppst,dec_ppst,roll_ppst


;Add the star tracker field - Spacecraft roll to the roll angles
;There is a 90degrees+17degres offset between spacecraft 
;coord and star tracker coord
  
  if n_params() lt 5 then begin
     print,'Please enter correct parameters: roll_check,year,doy,ra,dec,roll'
     return
  endif
  
  roll_start=roll_ppst
  roll_ppst=253.0-roll_ppst
  
;Read the star catalog
  readcol,'/home/pagani/XDS/star_track/test_shuffle/fl_cat.txt',ra_cat,dec_cat,mag_cat,id_cat,format='(f,f,f,a)',skip=1,/silent
  
;DENSITY
  radec6deg=0
  radecdens=0
  border=0

;Check for stars in FoV 
  
  for inccat=0l,n_elements(ra_cat)-1 do begin
     res = sphdist(ra_cat(inccat), dec_cat(inccat), ra_ppst, dec_ppst,/degrees) ;angles in DEGREES
     if res lt 6.0 then begin                                                   ; if star closer than 6.0 degrees check if it's in the 8x8 camera field of view
        radec6deg=radec6deg+1
        rotate_roll,ra_ppst,dec_ppst,ra_cat(inccat),dec_cat(inccat),roll_ppst*!dtor,ra_cat_new,dec_cat_new
                                ;Get Star Tracker H and V coordinates after rotation
        
                                ;Here I do a simple check if there is
                                ;a situation like this: RA_PPST+359.8
                                ;and RA_CAT=0.5, then the two sources
                                ;are very close but H=359!!! not real!!
        if (abs(ra_cat_new-ra_ppst) gt 340.) then begin
           if ra_cat_new gt ra_ppst then ra_cat_new=ra_cat_new-360. else ra_cat_new=ra_cat_new+360
        endif
        H=(ra_cat_new-ra_ppst)*cos(dec_ppst*!dtor) ;Star Tracker H,V stars coord
        V=dec_cat_new-dec_ppst
        
                                ;Checks if star in FoV and if it's at the border
        
        if abs(H) le 4.0 and abs(V) le 4.0 then begin
           prvar='H= '+strmid(strtrim(string(H),2),0,5)+'  V= '+strmid(strtrim(string(V),2),0,5)+'   Mag='+strmid(strtrim(string(mag_cat(inccat)),2),0,5)+'    ID: '+strtrim(string(id_cat(inccat)),2)
           print,prvar
           radecdens=radecdens+1 ; That is if the star is in the field of view
           if  abs(H) gt 3.75 or abs(V) gt 3.75 then border=border+1
        endif
        
                                ;if abs(H) ge 4.0 or abs(V) ge 4.0 then begin
                                ;    print,'!!!Star just outside field of view!!!   H= ',H,'  V=',V,'   Mag=',mag_cat(inccat),'    ID:',id_cat(inccat)
                                ;    print,ra_cat(inccat),dec_cat(inccat)
                                ;endif
        
    endif
  endfor
  
;If problems, print warning
  
  if (radecdens-border le 2) then print,'WARNING!!!  Stars at border causing problems!!! Not enough stars in FoV'
  print,'Number of stars within pointing 6 degrees / in star tracker Field of View:  ',radec6deg,radecdens
  
  
  ;CALCULATE ALLOWED ROLL ANGLE VALUES.  I have 2 options, one is to use Jamie's roll_range.py, the other is to calculate the optimal roll angle and get a range of +/-9.5 around the optimal roll 
  
  if roll_start lt 360. then begin

     print,''
     print,'Running roll_range.py to find allowed roll angles....'
     print,''
     rollcomando='roll_range.py '+strtrim(string(year),2)+' '+strtrim(string(doy),2)+' '+strtrim(string(ra_ppst),2)+' '+strtrim(string(dec_ppst),2)+' > roll_log.txt'
     print,rollcomando
;   stop
     spawn,rollcomando
     readcol,'roll_log.txt',roll_allowed,format='(f)',/silent
     spawn,'rm -f roll_log.txt'
   
                                ;Jamie's roll range starts from the optimal roll range and then it goes up or down according to temperature
                                ;if it goes up 15 degrees the allowed range will be optimal +/- 15 degrees.
                                ;Here I get the maximum and minimum roll angle, taking care of the 0/360 problem
   
     optimal=roll_allowed[0]
     endroll=roll_allowed[n_elements(roll_allowed)-2]
     if (endroll-optimal) gt 0. then begin
        maxroll=endroll
        minroll=optimal-(maxroll-optimal)
     endif else begin
        minroll=endroll
        maxroll=optimal+(optimal-minroll)
     endelse
     
     
  endif 
  
  if roll_start ge 360 then begin
     
     ;Get Sun RA,DEC
     
     ydn2md, year, doy, mese, giorno
     jdcnv, year, mese, giorno, 0 ,jd ;Find Julian date jd = 2445090.5   
     sunpos, jd, rasun, decsun
     
     ;Get Swift optimal roll_angle
     
     vector,decsun,rasun,vSun
     vector,dec_ppst,ra_ppst,vT
  ;  # get cross product of Sun vector and target vector (vSun x vT) :  result is vector Y
     vY = crossp(vT, vSun)
     vnY=vY/norm(vY)
   ; # get cross product of normalized Y and target vector: result is vector Z
     vZ=  crossp(vnY, vT)
     vnZ=vZ/norm(vZ)

     if (vnY[2] ne 0) and (vnZ[2] ne 0) then begin
        newroll = atan(vnY[2],vnZ[2])
        newroll = newroll/!dtor
     endif else begin
;        # roll is not uniquely defined, arbitrarily pick 0.0 or 180.0
        newroll = 0
        if (vSun[0]*(-cos(newroll)*sin(dec_ppst*!dtor)*cos(ra_ppst*!dtor) - sin(newroll)*sin(ra_ppst*!dtor)) + vSun[1]*(-cos(newroll)*sin(dec_ppst*!dtor)*sin(ra_ppst*!dtor) + sin(newroll)*cos(ra_ppst*!dtor)) + vSun[2]*cos(newroll)*cos(dec_ppst*!dtor) < 0.0) then newroll=180.0
                                ;          # 0.0 would put the sun vector in -Z, so use 180.0 deg
     endelse 
     if (newroll lt 0.0) then newroll = newroll+360.0
                                ;return newroll*!dtor
     
     minroll=newroll-9.5
     maxroll=newroll+9.5
  endif
  
          
  if minroll gt 360 then rollprintamin=minroll-360 else rollprintamin=minroll 
  if minroll lt 0 then rollprintamin=minroll+360 else rollprintamin=minroll 
  
  if maxroll gt 360 then rollprintamax=maxroll-360 else rollprintamax=maxroll 
  if maxroll lt 0 then rollprintamax=maxroll+360 else rollprintamax=maxroll 
  
  print,'Allowed roll angles:',rollprintamin,rollprintamax
  openw,report,'roll_report.txt',/get_lun
  printf,report,'Allowed roll angles:',rollprintamin,rollprintamax
  rolluse=253.0-minroll
  roll_stars=intarr(round(maxroll-minroll)+1)
  roll_array=fltarr(round(maxroll-minroll)+1)
  for cont_roll=0,round(maxroll-minroll) do begin
     radec6deg=0
     radecdens=0
     border=0
     for inccat=0l,n_elements(ra_cat)-1 do begin
        res = sphdist(ra_cat(inccat), dec_cat(inccat), ra_ppst, dec_ppst,/degrees)  ;angles in DEGREES
        if res lt 6.0 then begin                                                    ; if star closer than 6.0 degrees check if it's in the 8x8 camera field of view
           radec6deg=radec6deg+1
           rotate_roll,ra_ppst,dec_ppst,ra_cat(inccat),dec_cat(inccat),rolluse*!dtor,ra_cat_new,dec_cat_new
                                ;stop
                                ;Get Star Tracker H and V coordinates after rotation
            
                                ;Here I do a simple check if there is
                                ;a situation like this: RA_PPST+359.8
                                ;and RA_CAT=0.5, then the two sources
                                ;are very close but H=359!!! not real!!
           if (abs(ra_cat_new-ra_ppst) gt 340.) then begin
               if ra_cat_new gt ra_ppst then ra_cat_new=ra_cat_new-360. else ra_cat_new=ra_cat_new+360
            endif
            H=(ra_cat_new-ra_ppst)*cos(dec_ppst*!dtor) ;Star Tracker H,V stars coord
            V=dec_cat_new-dec_ppst
            
        ;Checks if star in FoV and if it's at the border
            
            if abs(H) le 4.0 and abs(V) le 4.0 then begin
               prvar='H= '+strmid(strtrim(string(H),2),0,5)+'  V= '+strmid(strtrim(string(V),2),0,5)+'   Mag='+strmid(strtrim(string(mag_cat(inccat)),2),0,5)+'    ID: '+strtrim(string(id_cat(inccat)),2)
               printf,report,prvar
               radecdens=radecdens+1 ; That is if the star is in the field of view
               if  abs(H) gt 3.75 or abs(V) gt 3.75 then border=border+1
            endif
            
                                ;if abs(H) ge 4.0 or abs(V) ge 4.0 then begin
        ;    print,'!!!Star just outside field of view!!!   H= ',H,'  V=',V,'   Mag=',mag_cat(inccat),'    ID:',id_cat(inccat)
        ;    print,ra_cat(inccat),dec_cat(inccat)
        ;endif

         endif
     endfor
     if minroll gt 360 then rollprinta=minroll-360 else rollprinta=minroll 
     if minroll lt 0 then rollprinta=minroll+360 else rollprinta=minroll  
     printf,report,'Roll: ',+strtrim(string(rollprinta),2)+' Stars in FoV : ',strtrim(string(radecdens-border),2) 
     rolluse=rolluse-1.
     roll_array(cont_roll)=minroll
     minroll=253.0-rolluse
     roll_stars(cont_roll)=radecdens-border
  endfor 
  
  maxstars=max(roll_stars)
  indexma=where(roll_stars eq maxstars)
  rangemin=min(roll_array(indexma))
  rangemax=max(roll_array(indexma))
  
  if rangemin gt 360 then rangemin=rangemin-360
  if rangemax gt 360 then rangemax=rangemax-360
  
  if rangemin lt 0 then rangemin=rangemin+360
  if rangemax lt 0 then rangemax=rangemax+360
  
  print,'Best Roll Angle range: '+strtrim(string(rangemin),2)+' - ',strtrim(string(rangemax),2)+' with ',strtrim(string(maxstars),2)+' stars.'
  printf,report,'Best Roll Angle range: '+strtrim(string(rangemin),2)+' - ',strtrim(string(rangemax),2)+' with ',strtrim(string(maxstars),2)+' stars.'
  free_lun,report
  spawn,'emacs roll_report.txt&'
  
end

