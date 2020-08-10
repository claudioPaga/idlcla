function Radiator_Earth_Angle, swift_p, swift_q
;Returns angle of radiation vs earth given a swift position vector and
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
  return, acos(-vbn[2])
end

pro plot_earth_angle_maxbrte, table_file

;Makes 2 histograms: 
;1) Radiator angle vs earth center at the min elv angle while in sunlight
;2) Radiator angle vs earth center at the peak of the bright earth cr
DEVICE,DECOMPOSED=0.
tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1

readcol, table_file, obsid, orbit, orb_start, orb_end, roll_angle, sun_angle, min_elv, min_elv_with_sun, min_elv_with_sun_time, min_elv_brte_count, max_brte_count, time_of_max_brte, elv_at_max_brte, sun_when_max, long_at_start, lat_at_start, long_at_end, lat_at_end, format='(a,i,d,d,i,i,i,i,d,l,l,d,i,i,i,i,i,i)'

;Here I derive the earth angles at time of cr peak
earth_angle_array_minelv38 = fltarr(n_elements(obsid))
earth_angle_array_max_brte = fltarr(n_elements(obsid))

openw,lunminelv,'table_swift_sao_earthangle_minelv.txt',/get_lun
openw,lunpeakcr,'table_swift_sao_earthangle_peakcr.txt',/get_lun

for i = 0, n_elements(obsid)-1 do begin
   name_saofile = '/swift/data/archive/obs/'+obsid[i]+'/auxil/sw'+obsid[i]+'sao.fits.gz'
   quat_file_match =  file_search(name_saofile, count =count_coord)
;   if count_coord eq 1 then begin ; Derive values if sao file with position and quaternion is found
   quat_file = mrdfits(name_saofile ,1,/silent)
   ;*************************************************
                                ;Finds the Swift pos/quat at time of
                                ;peak or count rate during an orbit
   index_minelv = min(abs(quat_file.TIME - min_elv_with_sun_time[i]), minelv_index)
   swift_p = quat_file(minelv_index).position
   swift_q = quat_file(minelv_index).quaternion
   earth_angle =  Radiator_Earth_Angle(swift_p, swift_q)
   earth_angle_array_minelv38[i] = earth_angle*180./!PI
                                ;  endif else earth_angle_array = 1000. ;If we filed to find quat/position file then just give 1000. to the angle
   printf,lunminelv, obsid[i], orbit[i], swift_p[0], swift_p[1], swift_p[2], swift_q[0], swift_q[1], swift_q[2], swift_q[3],earth_angle_array_minelv38[i], format='(a13,i3,d13.5,d13.5,d13.5,d12.7,d12.7,d12.7,d12.7,d8.2)'
      
                                ;*************************************************
                                ;Finds the Swift pos/quat at time of
                                ;min elv angle when sunlit
   index_brte_peak = min(abs(quat_file.TIME - time_of_max_brte[i]), minpeak_index)
   swift_p = quat_file(minpeak_index).position
   swift_q = quat_file(minpeak_index).quaternion
   earth_angle =  Radiator_Earth_Angle(swift_p, swift_q)
   earth_angle_array_max_brte[i] = earth_angle*180./!PI
 ;  endif else earth_angle_array = 1000. ;If we filed to find quat/position file then just give 1000. to the angle
   printf,lunpeakcr, obsid[i], orbit[i], swift_p[0], swift_p[1], swift_p[2], swift_q[0], swift_q[1], swift_q[2], swift_q[3],earth_angle_array_max_brte[i], format='(a13,i3,d13.5,d13.5,d13.5,d12.7,d12.7,d12.7,d12.7,d8.2)'
endfor

free_lun,lunminelv
free_lun,lunpeakcr

brte_cutoff = 50
;*****************************************
;Select orbits with elv<38 (and sunlit), and orbits with elv<38 that
;are bright at elv<38
elv38_index = where(min_elv_with_sun le 38)
elv38_earth_angle = earth_angle_array_minelv38[elv38_index] ;plothist 1, 1st array
elv38_br_index = where(min_elv_with_sun le 38 and min_elv_brte_count ge brte_cutoff)
elv38_br_earth_angle = earth_angle_array_minelv38[elv38_br_index] ;plothist 1, 2nd array
;*********************************************
;Select orbits with elv<38 (and sunlit), and orbits with elv<38 that
;are bright at ANY elv angle 
peak_earth_angle = earth_angle_array_max_brte[elv38_index]
peak_br_index = where(min_elv_with_sun le 38 and  max_brte_count ge brte_cutoff)
peak_br_earth_angle = earth_angle_array_max_brte[peak_br_index]

;PLOTS
;**************************
!p.multi= [0,2,2,1,1]
nameps= 'earth_angle_hist.ps'
set_plot,'PS'
device,filename=nameps,    /color,/landscape ;, xs=30,ys=20

bin_variable = 10.
;Plot histo of earth_angle for orbits with ELV<38 and and all bright
;orbits with ELV<38 which are bright at ELV<38
earth_38_histo = histogram(elv38_earth_angle, bin = bin_variable, min = 0, locations = xearthhisto_38)
earth_brte_38_histo = histogram(elv38_br_earth_angle, bin = bin_variable, min = 0, locations = xearthhisto_38_brte)
;plothist, elv38_earth_angle, xearth38, yearth38, background = 1, color = 0, ytit='n orbits', xtit= 'earth angle', tit = 'earth angle at min(ELV) in sunlight', bin = bin_variable, xr = [0,180],xstyle = 1, yr = [0,max_yearth], charsize = 1.1
;max_yearth = max(yearth38)
plothist, elv38_earth_angle, xearth38, yearth38, background = 1, color = 0, ytit='n orbits', xtit= 'earth angle', tit = 'earth angle at min(ELV) in sunlight', bin = bin_variable, xr = [0,180],xstyle = 1, charsize = 1.05
plothist, elv38_br_earth_angle, xearth_brte38, yearth_brte38, /overplot, color = 3, /fill, fcolor = 3, bin = bin_variable
;Plot of the ratio
plot, xearthhisto_38, earth_brte_38_histo*1./earth_38_histo, background = 1, color = 0, psym=2, tit = 'fraction of orbits affected by bright earth', xtit = 'earth angle', ytit = 'fraction affected by bright earth',xstyle = 1, charsize = 1.05, xr = [0,180], yr = [0,1]

;Plot histo of earth_angle for orbits with ELV<38 and and all bright
;orbits with ELV<38 which are bright at ANY ELV
;brte_elv38_earth_angle_fullorb = [0, brte_elv38_earth_angle_fullorb]
earth_38_histo_fullorb = histogram(peak_earth_angle, bin = bin_variable, min = 0, locations = xearthhisto_38_fullorb)
earth_brte_38_histo_fullorb = histogram(peak_br_earth_angle, bin = bin_variable, min = 0, locations = xearthhisto_38_brte_fullorb)
;max_yearth_full = max(earth_38_histo_fullorb)
plothist, peak_earth_angle, xearth38_fullorb, yearth38_fullorb, background = 1, color = 0, ytit='n orbits', xtit= 'earth angle', bin = bin_variable, xr = [0,180],xstyle = 1, charsize = 1.05, tit = 'earth angle at bright earth peak'
plothist, peak_br_earth_angle, xearth_brte38_fullorb, yearth_brte38_fullorb, /overplot, color = 3, /fill, fcolor = 3, bin = bin_variable
;Plot of the ratio
plot, xearthhisto_38_fullorb, earth_brte_38_histo_fullorb*1./earth_38_histo_fullorb, background = 1, color = 0, psym=2, xtit = 'earth angle', ytit = 'fraction affected by bright earth',xstyle = 1, xr = [0,180], yr = [0,1], charsize = 1.05, tit = 'fraction of orbits affected by bright earth'

;screen_text1 = 'all orbits'
;XYOUTS, 0.07, 0.95, screen_text1, /norm, color = 0, charsize = 0.8
;screen_text2 = 'all orbits with bright earth(BE)'
;XYOUTS, 0.07, 0.93, screen_text2, /norm, color = 3, charsize = 0.8
;screen_text3 = 'orbits in eclipse'
;XYOUTS, 0.4, 0.835, screen_text3, /norm, color = 2, charsize = 0.8
screen_text1 = 'sunlit orbits (elv<38)'
XYOUTS, 0.15, 0.93, screen_text1, /norm, color = 0, charsize = 0.8
screen_text2 = 'sunlit orbits with BE @ elv<38'
XYOUTS, 0.15, 0.90, screen_text2, /norm, color = 3, charsize = 0.8
screen_text1 = 'sunlit orbits (elv<38)'
XYOUTS, 0.65, 0.93, screen_text1, /norm, color = 0, charsize = 0.8
screen_text2 = 'sunlit orbits (elv<38) with BE'
XYOUTS, 0.65, 0.90, screen_text2, /norm, color = 3, charsize = 0.8

device, /close
set_plot,'x'
;WRITE_PNG, 'earth_angle_brte_hist_plot_with_elv38.png', /transparent, tvrd(true =1)

end
