function Radiator_Earth_Angle, swift_p, swift_q

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


pro derive_earth_angle, table_file
;
;CP, 20/Aug/2012
;
;Reads in info from an orbit, and from the Quaternion and the Swift
;coordinates (in Km) derives the angle of Swift with respect to Earth
;

;

if N_PARAMS() ne 1 then begin
    print, 'WARNING!!!!!!'
    print, 'Please enter the obsids table file'
    print, "Ex: IDL> derive_earth_angle, 'grbs_obsids.txt'"
    print, '!!!!!!!!!!!!!'
    return
endif

DEVICE,DECOMPOSED=0.
tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1

readcol, table_file, obsid, orbit, orb_start, orb_end, roll_angle, sun_angle, min_elv, min_elv_with_sun, min_elv_with_sun_time, min_elv_brte_count, max_brte_count, time_of_max_brte, elv_at_max_brte, sun_when_max, long_at_start, lat_at_start, long_at_end, lat_at_end, format='(a,i,d,d,i,i,i,i,d,l,l,d,i,i,i,i,i,i)'

for i = 0, n_elements(obsid)-1 do begin
   name_saofile = '/swift/data/archive/obs/'+obsid[i]+'/auxil/sw'+obsid[i]+'sao.fits.gz'
   quat_file_match =  file_search(name_saofile, count =count_coord)
   if count_coord eq 1 then begin ; Derive values if sao file with position and quaternion is found
      quat_file = mrdfits(name_saofile ,1,/silent)

   ;Finds the Swift pos/quat at time of br. earth peak
      index_brte_peak = min(abs(quat_file.TIME - time_of_max_brte), min_index)
      
      swift_p = quat_file(min_index).position
      swift_q = quat_file(min_index).quaternion
      swift_orig =  swift_q

      swift_q[0] =  swift_orig[3]
      swift_q[1:3] =  swift_orig[0:2]
      earth_angle =  Radiator_Earth_Angle(swift_p, swift_orig)
      earth_angle_deg = earth_angle*180./!PI
      print,  'Earth angle with original quaternion', earth_angle_deg
      
      earth_angle =  Radiator_Earth_Angle(swift_p, swift_q)
      earth_angle_deg = earth_angle*180./!PI
      print,  'Earth angle with modified quaternion', earth_angle_deg
      
   endif
endfor
end
