function Radiator_Earth_Angle, swift_p, swift_q
;Returns angle of radiator vs earth given a swift position vector and
;a quaternion vector

 ; swift_q is the spacecraft quaternion. If q = (0, 0, 0, 1) (ie no rotation)
 ; then the spacecraft X axis (XRT pointing direction) points to the
 ; vernal equinox and the Z axis to the celestial north pole (or the
 ; radiator normal, -Z, points to celestial south).
 ; va is the negative of the spacecraft position vector in Earth coords.
 ; vb is the earth position vector in spacecraft coords  

  va = -swift_p
  q1 = swift_q[0]
  q2 = swift_q[1]
  q3 = swift_q[2]
  q4 = swift_q[3]
  
  qmag = sqrt(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4)
  if (qmag > 0.0001) then begin
     q1 = q1 / qmag
     q2 = q2 / qmag
     q3 = q3 / qmag
     q4 = q4 / qmag
  endif else begin
     print, "*** Warning: qmag < 0.0001 ***"
     q1 = 0.0
     q2 = 0.0
     q3 = 0.0
     q4 = 1.0
  endelse

  w1 = q2 * va[2] - q3 * va[1]
  w2 = q3 * va[0] - q1 * va[2]
  w3 = q1 * va[1] - q2 * va[0]

  u1 = q2 * w3 - q3 * w2
  u2 = q3 * w1 - q1 * w3
  u3 = q1 * w2 - q2 * w1

  vb1 = va[0] + 2 * (u1 - q4 * w1)
  vb2 = va[1] + 2 * (u2 - q4 * w2)
  vb3 = va[2] + 2 * (u3 - q4 * w3)

  vb = [vb1, vb2, vb3]
  vbmag = sqrt(vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2])
  vbn = [vb[0] / vbmag , vb[1] / vbmag , vb[2] / vbmag]
                                ; # radiator normal is on the spacecraft -Z axis
 ;   # angle of radiator normal to earth is the dot product of
 ;   # vbn with (0, 0, -1)
 ;   # which is just -vbn[2]
  return, acos(-vbn[2])*180./!PI
end
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pro gti_check, xrt_gti, xrt_gti_new

  gti_new_start = xrt_gti.start
  gti_new_stop = xrt_gti.stop
  if n_elements(gti_new_start) gt 1 then begin

;clearing out fake orbits
     for i = 1, n_elements(gti_new_start)-1 do begin
re_clear:
        if  gti_new_start[i] ne 0 then begin
           if gti_new_start[i] - gti_new_stop[i-1] lt 2000. then begin 
              for j = i, n_elements(gti_new_start)-1 do begin
                 gti_new_stop[j-1] = gti_new_stop[j]
              endfor
              for  j = i, n_elements(gti_new_start)-2 do begin
                 gti_new_start[j] = gti_new_start[j+1]
              endfor
              gti_new_start[n_elements(gti_new_start)-1]=0
              gti_new_stop[n_elements(gti_new_stop)-1]=0
              goto, re_clear
           endif
        endif
     endfor 
     
     gti_new_start_no0_index= where(gti_new_start ne 0, n_no0)
; Create final structure with cleaned GTIS(no fake orbits)
     xrt_gti_new = {start: 0d, stop: 0d}
     xrt_gti_new = replicate(xrt_gti_new,n_no0)
     xrt_gti_new.start = gti_new_start[gti_new_start_no0_index]
     xrt_gti_new.stop = gti_new_stop[gti_new_start_no0_index]
  endif else xrt_gti_new = xrt_gti
end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro plot_radiator_earth_angle_satfile, filename_evt, saofile, mkffile

;CP, 10 May 2013
;
;PURPOSE: check for possible bright earth contamination based on the
;elevation angle and the radiator earth angle
;We now have this constraint on board: IF ELV<40 ==> RAD-EARTH >50
;But sometimes this ELV<40 is not enough, we see some contamination
;even at ELV between 40 and 50, so the best condition should have been
;ELV<50.
;
;So, in these plots, if ELV<50 and RAD-EARTH<50 then BRIGHT EARTH is possible!! 
;
;Given a swift event file, sao.fits file and mkf fileplot the radiator vs earth angle during
;the observations
;
;

;Ex:plot_radiator_earth_angle_satfile,'sw00050400128xpcw3po_uf.evt','sw00050400128sao.fits','sw00050400128s.mkf'

tab = mrdfits(saofile,1,hd1)
swift_mkf = mrdfits(mkffile,1,/silent)
xrt_gti = mrdfits(filename_evt,2,/silent)
                                ;REMOVES SHORT GAPS IN GTIS
gti_check, xrt_gti, xrt_gti_check

time = tab.TIME
swift_p = tab.position
swift_q = tab.quaternion

earth_angle_deg = fltarr(n_elements(time))

for j = 0l, n_elements(time)-1 do earth_angle_deg[j] =  Radiator_Earth_Angle(swift_p[*,j], swift_q[*,j])

index = where(Finite(earth_angle_deg) ne 0)

obsid = strmid(filename_evt,3,11)

DEVICE,DECOMPOSED=0.
tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1

nameps= 'earthangle_brightearth_nocr'+strtrim(obsid,2)+'.ps'
set_plot,'PS'
;opens file
device,filename=nameps,    /color,/landscape ;, xs=30,ys=20

for k = 0, n_elements(xrt_gti_check.start)-1 do begin
    minplot =300
   maxplot = 0

   indexelv = where(swift_mkf.TIME gt xrt_gti_check[k].start and swift_mkf.TIME lt xrt_gti_check[k].stop)
   maxelv = max(swift_mkf[indexelv].ELV)
   minelv = min(swift_mkf[indexelv].ELV)
   if minelv lt minplot then minplot = minelv
   if maxelv gt maxplot then maxplot = maxelv

   indexearth = where(tab.TIME gt xrt_gti_check[k].start and tab.TIME lt xrt_gti_check[k].stop,nindexearth)
   maxearth = max(earth_angle_deg[indexearth])
   minearth = min(earth_angle_deg[indexearth])
   if minearth lt minplot then minplot = minearth
   if maxearth gt maxplot then maxplot = maxearth
   
   if minelv lt 50 and minearth lt 50 then be_condition = 1 else be_condition = 0

   titolo = 'Orbit '+strtrim(string(k+1),2)
   plot,swift_mkf[indexelv].TIME, swift_mkf[indexelv].ELV , psym=3, color=0, background=1, xr=[xrt_gti_check[k].start,xrt_gti_check[k].stop], yr=[minplot,maxplot],tit=titolo
   xyouts,0.5,0.9,'Possible BE if ELV<50 AND RAD-EARTH <50',/norm,color=0
   if be_condition then xyouts,0.5,0.85,'WARNING!!!!, POSSIBLE BE CONTAMINATION IN THIS ORBIT!!!!',/norm,color=3
   xyouts,0.8,0.8,'ELV',/norm,color=0
   if nindexearth ge 2 then begin
       oplot, tab[indexearth].TIME, earth_angle_deg[indexearth], psym=3, color=2
       xyouts,0.8,0.75,'RADIATOR-EARTH',/norm,color=2
       ;indexlow = where(earth_angle_deg[indexearth] lt 60., numero)
       ;if numero gt 1 then begin
       ;    t1 = tab[indexearth].TIME
       ;    print, t1[indexlow[0]], format='(d15.1)'
       ;endif
   endif
   xlimit = [xrt_gti_check[k].start-100000.,xrt_gti_check[k].stop+100000.]
   ylimit = [50,50]
   oplot, xlimit, ylimit, color=3
endfor


 device,/close
 spawn,'gv '+nameps+'&'
  set_plot,'x'
  
end
