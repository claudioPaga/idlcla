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
  return, acos(-vbn[2])
end


pro plotting_routine, obsid, xrt_gti, time_br, raobj, decobj, papnt, object, swift_mkf, xrt_dataset_unfilt, swift_coordinates, swift_sao, earth_angle

  nameps= 'earthangle_brightearth'+strtrim(string(obsid),2)+'.ps'
  set_plot,'PS'
;opens file
  device,filename=nameps,    /color,/landscape ;, xs=30,ys=20

  hist_bin = 20

;~~~~~~~~~~~~~
;variables for counting brte orbs
  number_of_brte_orbs = 0
  number_of_brte_orbs_38_sun = 0
;~~~~~~~~~~~~~~
;array for finding out how many orbits go below 38 elv
  br_orb_min_elv_array = fltarr(n_elements(xrt_gti.start)) + 360

;Selection to exclude very short orbits
  if n_elements(xrt_gti.stop) gt 1 then begin
     min_gti_lenght = 50.
     gti_lenght = xrt_gti.stop - xrt_gti.start
     lenght_index = where(gti_lenght gt min_gti_lenght)
     gti_long_start = xrt_gti[lenght_index].start
     gti_long_stop = xrt_gti[lenght_index].stop
  endif else begin
     gti_long_start = xrt_gti.start
     gti_long_stop = xrt_gti.stop
  endelse
     

;PRCESS DATA FOR EACH ORBIT

;ORBITS LOOP START
  for j= 0,n_elements(gti_long_start)-1 do begin

;SELECT HK FROM MKF FILE FOR EACH ORBIT
     orbit_index = where(swift_mkf.time ge gti_long_start[j] and swift_mkf.time le gti_long_stop[j])
     orbit= swift_mkf[orbit_index]
     orbit_sun_index = where(orbit.sunshine eq 1, n_sunshine) 
     if n_sunshine gt 0 then orbit_sun = orbit[orbit_sun_index]
     orbit_rollangle = median(orbit.roll)
     orbit_sunangle = median(orbit.sun_angle)
     swift_mkf_elv_min = min(orbit.elv)
     
;PLOTTING START
     m= j mod 8
     if m eq 0 then !p.multi= [0,4,2] else !p.multi= [8-m,4,2]
     histogram_time_br= histogram(time_br, bin=hist_bin)
     if n_elements(histogram_time_br) gt 1 then begin ;Makes histo if there are enough databins, otherwise just plot zeros
        plothist, time_br-gti_long_start[j], xr=[0,gti_long_stop[j]-gti_long_start[j]],/ylog, yr=[0.1, max(histogram_time_br)], bin=20,xtit='time', ytit='counts', tit= '', background= 1, color=0, ystyle=8
     endif else begin
        xpl = indgen(fix(gti_long_stop[j]-gti_long_start[j]))
        ypl = intarr(n_elements(xpl))
        plot, xpl, ypl, xr=[0,gti_long_stop[j]-gti_long_start[j]],/ylog, yr=[0.1, 10.], xtit='time', ytit='counts', tit= '', background= 1, color=0, ystyle=8
     endelse

                           ;start of displaying info about the OBSID
     screen_text1 = 'ObsID = '+strtrim(string(obsid),2)
     XYOUTS, 0.02, 0.49, screen_text1, /norm
     screen_text = '(RA, DEC) = '+strtrim(string(raobj),2)+','+strtrim(string(decobj),2)
     XYOUTS, 0.22, 0.49, screen_text, /norm
     screen_text2 = 'roll angle = '+strtrim(string(papnt),2)
     XYOUTS, 0.5, 0.49, screen_text2, /norm
     screen_text3 = 'object name = '+strtrim(string(object),2)
     XYOUTS, 0.7, 0.49, screen_text3, /norm
                                ;end
     
;Overplot the ELV angle
     !p.multi= [8-m,4,2]  
     plot, swift_mkf.time-gti_long_start[j], swift_mkf.elv, color=0, line=2, /noerase, yrange=[0,180], ystyle=4, xrange=!x.crange
     axis, yaxis=1, yrange=!y.crange, ytit='_ _ _ _ _elv angle', color=0, xstyle=-1
;Overplot the ELV angle in red during SUNSHINE times
     swift_mkf_sun_index = where(swift_mkf.sunshine eq 1, n_sunshine_plot) 
     if n_sunshine_plot gt 1 then oplot, swift_mkf[swift_mkf_sun_index].time-gti_long_start[j], swift_mkf[swift_mkf_sun_index].elv, color=3, psym=3
                                ;update indices to select plotting quadrant
     k = j mod 4
     h= j mod 8
     if h le 3  then f = 1 else f = 0
    
;Insert plot of Swift coordinates during the orbit
;     coord_orbit_index = where(swift_coordinates.time ge gti_long_start[j] and swift_coordinates.time le gti_long_stop[j])
;     if  n_elements(swift_coordinates[coord_orbit_index].position[0]) gt 1 then plot, swift_coordinates[coord_orbit_index].position[0], swift_coordinates[coord_orbit_index].position[1], position = [(0.25*k)+0.075,0.35+(f*0.5),(0.25*k)+0.23,0.45+(f*0.5)], /noerase, color = 4, xr= [0,360], yr = [-25,25], xstyle = 1, ystyle = 1, psym = 3

;Insert plot of Swift earth angle during the orbit
;     earth_orbit_index = where(swift_sao.time ge gti_long_start[j] and swift_sao.time le gti_long_stop[j])
;     if  n_elements(swift_sao[earth_orbit_index].time) gt 1 then plot, swift_sao[earth_orbit_index].time-gti_long_start[j], earth_angle[earth_orbit_index], position = [(0.25*k)+0.075,0.35+(f*0.5),(0.25*k)+0.23,0.45+(f*0.5)], /noerase, color = 4, xr=[0,gti_long_stop[j]-gti_long_start[j]], yr = [0,180], xstyle = 1, ystyle = 1, psym = 3

;Overplot the Radiator normal angle
     earth_orbit_index = where(swift_sao.time ge gti_long_start[j] and swift_sao.time le gti_long_stop[j], nindex_earth)
     if  nindex_earth gt 1 then begin
        plot, swift_sao[earth_orbit_index].time-gti_long_start[j], earth_angle[earth_orbit_index], color=0, /noerase, yrange=[0,180], ystyle=4, xrange=!x.crange
        oplot, swift_sao[earth_orbit_index].time-gti_long_start[j], earth_angle[earth_orbit_index], color=4, psym=3
        axis, yaxis=1, yrange=!y.crange, ytit='_ _ _ _ _elv angle', color=0, xstyle=-1
     endif

;PLOTTING END
     
;here I print the start and stop time when there should not be be, ie,
;good gtis to use to remove be (mostly for pls)

;index_gti_elv_pls = where(swift_mkf.elv gt 40., n_elv_pls)
;index_gti_eartj_pls = where(earth_angle[earth_orbit_index] gt 50., n_earth_pls)
;if n_elv_pls gt 0 and n_earth_pls gt 0 then begin
 ;  time_elv_ok_pls = swift_mkf[index_gti_elv_pls].time
 ;  time_earth_check = swift_sao[earth_orbit_index].time
 ;  time_earth_ok_pls = time_earth_check[index_gti_eartj_pls]
  ; pls_tstart_evl = min(time_elv_ok_pls)
  ; pls_tstop_evl = max(time_elv_ok_pls)
  ; pls_tstart_earth = min(time_earth_ok_pls)
  ; pls_tstop_earth = max(time_earth_ok_pls)
;   timeends = [pls_tstop_evl,pls_tstop_earth]
 ;  print, gti_long_start[j], min(timeends)
;endif else print, 'No good times in Orbit ',j



;DERIVE RELEVANT INFORMATION OF EACH ORBIT
                
       ;DERIVE APPROPRIATE START_TIMES AND STOP_TIMES DURING ORBIT
       time_br_orb_index = where(time_br ge gti_long_start[j] and time_br le gti_long_stop[j], n_time_br_orb_index)
       time_br_unfilt_orb_index = where(xrt_dataset_unfilt.time ge gti_long_start[j] and xrt_dataset_unfilt.time le gti_long_stop[j])
       min_histo_time = xrt_dataset_unfilt[time_br_unfilt_orb_index[0]].time ;Time of first event in the orbit
       max_histo_time = xrt_dataset_unfilt[time_br_unfilt_orb_index[n_elements(time_br_unfilt_orb_index)-1]].time;Time of last event in the orbit
        
       ;DERIVE INFORMATION ABOUT THE ORBIT IF THE ORBIT IS NOT EMPTY
       if n_time_br_orb_index gt 0 then begin
                                ;selecting times in orbit
            time_br_orb = time_br[time_br_orb_index]
            histogram_time_br_orb= histogram( time_br_orb, bin=hist_bin, locations = time_of_bin_start, min=min_histo_time, max=max_histo_time)
                                ;matching the events within our orbits times
            orbit_index = where(swift_mkf.time ge time_of_bin_start[0] and swift_mkf.time le time_of_bin_start[n_elements(time_of_bin_start)-1]+hist_bin)
            orbit= swift_mkf[orbit_index]
                                ;identifying where the events with the sun hitting the satellite are
            orbit_sun_index = where(orbit.sunshine eq 1, n_sunshine) 
      
            ;DERIVE MAX COUNTRATE AND ITS TIME AND
            ;THE CORRESPONDING ELV, LONG and LAT DURING THE ORBIT
            histogram_time_br_orb_max = max(histogram_time_br_orb, indexmax_hist)
            time_of_max = time_of_bin_start[indexmax_hist]
            
            elv_at_max_value  = min(abs(orbit.time - time_of_max), elv_at_max_index)
            elv_at_max = orbit[elv_at_max_index].elv
            diff_time = abs(swift_coordinates.time - time_of_max)
            mintimediff = min(diff_time , longlat_index) 
            long_at_brte_max = swift_coordinates[longlat_index].position[0]
            lat_at_brte_max = swift_coordinates[longlat_index].position[1]
         
            if histogram_time_br_orb_max ge 20. then number_of_brte_orbs = number_of_brte_orbs+1 ;CONSIDER THE ORBIT AFFECTED BY BRIGHT EARTH IF PEAK COUNTRATE IS ABOVE 20
            
            ;CHECK IF IN SUNSHINE DURING PEAK COUNTRATE OF THE ORBIT
            orbit_when_max_index = where(orbit.time ge time_of_max and orbit.time le time_of_max + 10)
            orbit_when_max = orbit[orbit_when_max_index]
            sun_when_max_count = orbit_when_max[0].sunshine ;

            ;DERIVE TIME AND VALUE OF MINIMUM ELV
            ;DURING SUNSHINE AND THE CORRESPONDING COUNTRATE
            if n_sunshine gt 0 then begin
               orbit_sun = orbit[orbit_sun_index]
               swift_mkf_elv_min_with_sun = min(orbit_sun.elv)
               br_orb_min_elv_array[j]= swift_mkf_elv_min_with_sun ; Populate the array with the value of the min(ELV) IN SUNSHINE during orbit j
                                ;finding where the (with sun) minnimum elv is in terms of the structure
               swift_mkf_elv_min_with_sun_where = where(orbit_sun.elv eq swift_mkf_elv_min_with_sun)
                                ;eliminating double readings of the lowest elv value
               if  n_elements(swift_mkf_elv_min_with_sun_where) gt 0 then swift_mkf_elv_min_with_sun_where_1=swift_mkf_elv_min_with_sun_where[0] else swift_mkf_elv_min_with_sun_where_1=swift_mkf_elv_min_with_sun_where
               swift_mkf_elv_min_with_sun_time= orbit_sun[swift_mkf_elv_min_with_sun_where_1].time
               hist_orb_min_angle_time_index = where(time_of_bin_start le swift_mkf_elv_min_with_sun_time and time_of_bin_start ge swift_mkf_elv_min_with_sun_time - hist_bin )
               histogram_orb_min_angle_count_value = histogram_time_br_orb[hist_orb_min_angle_time_index]
               if  histogram_orb_min_angle_count_value ge 20. and swift_mkf_elv_min_with_sun lt 38. then number_of_brte_orbs_38_sun = number_of_brte_orbs_38_sun+1 ;COUNTS ORBITS WITH BRIGHT EARTH WHILE ELV<38.
            endif

            ;POPULATE RELEVANT VARIABLES AND
            ;ARRAYS IN CASE OF ECLIPSE DURING ORBIT
            if n_sunshine  eq 0 then begin
               swift_mkf_elv_min_with_sun = 360
               swift_mkf_elv_min_with_sun_time = 0
               histogram_orb_min_angle_count_value = 0.
            endif
            
         endif else begin ;END OF DERIVING ORBIT INFO IF ORBIT IS NOT EMPTY
            
            ;POPULATE ARRAYS IN CASE OF EMPTY ORBIT
            orbit_index = where(swift_mkf.time ge gti_long_start[j] and swift_mkf.time le gti_long_stop[j])
            orbit= swift_mkf[orbit_index]
            orbit_sun_index = where(orbit.sunshine eq 1, n_sunshine)
            if n_sunshine gt 0 then begin
               orbit_sun = orbit[orbit_sun_index]
               swift_mkf_elv_min_with_sun = min(orbit_sun.elv)
               br_orb_min_elv_array[j]= swift_mkf_elv_min_with_sun
               swift_mkf_elv_min_with_sun_where = where(orbit_sun.elv eq swift_mkf_elv_min_with_sun)
               if  n_elements(swift_mkf_elv_min_with_sun_where) gt 0 then swift_mkf_elv_min_with_sun_where_1=swift_mkf_elv_min_with_sun_where[0] else swift_mkf_elv_min_with_sun_where_1=swift_mkf_elv_min_with_sun_where
               swift_mkf_elv_min_with_sun_time= orbit_sun[swift_mkf_elv_min_with_sun_where_1].time
            endif
            
            ;SETTING VARIABLES TO DEFAULT VALUES FOR EMPTY ORBITS
            histogram_time_br_orb_max = 0.
            time_of_max = 0.
            elv_at_max = 0.
            histogram_orb_min_angle_count_value = 0.
            sun_when_max_count = 0.
            long_at_brte_max = 5000
            lat_at_brte_max = 5000
            if n_sunshine  eq 0 then begin
               swift_mkf_elv_min_with_sun = 360.
               swift_mkf_elv_min_with_sun_time = 0.
               histogram_orb_min_angle_count_value = 0.
            endif
         endelse ;END OF POPULATING ARRAYS DURING EMPTY ORBITS

         ;DERIVE SWIFT COORDINATES INFO DURING ORBIT
         if  n_elements(coord_orbit_index) gt 1 then begin
            long_start = swift_coordinates[coord_orbit_index[0]].position[0]
            lat_start = swift_coordinates[coord_orbit_index[0]].position[1]
            long_end = swift_coordinates[coord_orbit_index[n_elements(coord_orbit_index)-1]].position[0]
            lat_end = swift_coordinates[coord_orbit_index[n_elements(coord_orbit_index)-1]].position[1]
         endif else begin
            long_start = 1000
            lat_start = 1000
            long_end = 1000
            lat_end = 1000
         endelse
         
                                ;PRINT OUTPUT FILE FOR EACH ORBIT

         namedata= 'obsid'+strtrim(string(obsid),2)+'_orbit'+strtrim(string(j+1),2)+'.txt'
         openw,lun,namedata,/get_lun
         printf, lun, obsid, j+1, gti_long_start[j], gti_long_stop[j], orbit_rollangle, orbit_sunangle, swift_mkf_elv_min, swift_mkf_elv_min_with_sun, swift_mkf_elv_min_with_sun_time, histogram_orb_min_angle_count_value, histogram_time_br_orb_max, time_of_max, elv_at_max[0], sun_when_max_count,long_start, lat_start, long_end, lat_end, long_at_brte_max, lat_at_brte_max, format = '(a12,i3,d13.1,d13.1,i5,i5,i4,i4,d13.1,i7,i7,d13.1,i4,i3,i5,i5,i5,i5,i5,i5)'
         Free_lun, lun
        
      endfor

;PRINT OUTPUT FILE TO SUMMARIZE THE ORBITS OF THE OBSID
;counting number of orbs with elv under 38 + sun
  orbits_le38_index = where(br_orb_min_elv_array le 38, number_of_orbs_below_38_elv_sun)
  number_of_orbits = n_elements(gti_long_start)
  namedata = 'obsid'+strtrim(string(obsid),2)+'allorbits_data.txt'
  openw,lun,namedata,/get_lun   
  printf, lun, obsid, number_of_orbits, number_of_brte_orbs, number_of_orbs_below_38_elv_sun, number_of_brte_orbs_38_sun
  Free_lun, lun
  
  device,/close
;spawn,'gv '+nameps+'&'
  set_plot,'x'

end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function exlude_hot_pixels, xrt_dataset_crop_2dhis, threshold_3sig, xrt_dataset_crop

;exlude hot pixels
  brte_index = where(xrt_dataset_crop_2dhis le threshold_3sig)
  brte_index_2d =  ARRAY_INDICES(xrt_dataset_crop_2dhis, brte_index)
  brte_index_2d[1,*] = brte_index_2d[1,*]+200
  xrt_dataset_string =  strtrim(string(xrt_dataset_crop.detx),2)+strtrim(string(xrt_dataset_crop.dety), 2)
  brte_index_string =  strtrim(string(brte_index_2d[0,*]), 2)+strtrim(string(brte_index_2d[1,*]), 2)
  brte_index_string1d =  strarr(n_elements(brte_index_string))
  for i =  0,  n_elements(brte_index_string)-1 do brte_index_string1d[i]=   brte_index_string[i]
  match2, xrt_dataset_string,  brte_index_string1d,  sub1,  sub2
  indexmatch =  where(sub1 ge 0)
  return, indexmatch            ;this is the array with indices of bright earth events, without hot pixels
end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function threshold_deriver, xrt_dataset_crop_2dhis
;this works out the threshold to discriminate between hot pixels and brte

xrt_dataset_crop_his_no0_index = where(xrt_dataset_crop_2dhis ge 0)
xrt_dataset_crop_his_no0 = xrt_dataset_crop_2dhis(xrt_dataset_crop_his_no0_index)
xrt_dataset_crop_his_no0_index_sorted= sort(xrt_dataset_crop_his_no0)
xrt_dataset_crop_his_no0_sorted=xrt_dataset_crop_his_no0[xrt_dataset_crop_his_no0_index_sorted]
xrt_dataset_crop_his_no0_sorted_66index = float(n_elements(xrt_dataset_crop_his_no0_sorted)*(2./3.))
xrt_dataset_crop_his_no0_sorted_66 = xrt_dataset_crop_his_no0_sorted[fix(xrt_dataset_crop_his_no0_sorted_66index)]
threshold_3sig= 3*xrt_dataset_crop_his_no0_sorted_66

return, threshold_3sig
end

;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro cropping_data, xrt_dataset, min_detx, max_detx, min_dety, max_dety, xrt_dataset_crop, xrt_dataset_crop_2dhis, stopper

;SELECTS BRIGHT EARTH AREA
;xrt_dataset_croped_index = array indices that match the selected area
;xrt_dataset_crop = array of the selected area

xrt_dataset_croped_index = where(xrt_dataset.detx gt min_detx and xrt_dataset.detx lt max_detx and xrt_dataset.dety gt min_dety and xrt_dataset.dety lt max_dety, total_matches)
stopper = 0
if total_matches lt 2 then stopper=1 else begin ;check for bright earth events to avoid empty arrays
   xrt_dataset_crop = xrt_dataset(xrt_dataset_croped_index)
   xrt_dataset_almostcrop_2dhis = HIST_2D(xrt_dataset_crop.detx, xrt_dataset_crop.dety, max1 = 79, max2 =399)
   xrt_dataset_crop_2dhis = xrt_dataset_almostcrop_2dhis[min_detx:max_detx-1, min_dety:max_dety-1]
endelse
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

pro brte_read32, filename_evt, filename_mkf, filename_coord, filename_sao, g32, obsid, raobj, decobj, papnt, object, xrt_dataset32, xrt_gti_check, swift_mkf, xrt_dataset, stopper, swift_coordinates, swift_sao, earth_angle

;1;putting eventlist into structure in idl
  xrt_dataset = mrdfits(filename_evt, 1, header1,/silent)
  stopper=0
  if g32 then xrt_dataset_index = where(xrt_dataset.grade eq 32, n_eve) else xrt_dataset_index = where(xrt_dataset.grade le 32, n_eve)
  
  if n_eve lt 3 then stopper = 1 else begin
     xrt_dataset32 = xrt_dataset[xrt_dataset_index]
;1.5;collecting data to be shown on the graph
     obsid= sxpar(header1, 'OBS_ID')
     raobj= sxpar(header1, 'RA_OBJ')
     decobj= sxpar(header1, 'DEC_OBJ')
     papnt= sxpar(header1, 'PA_PNT')
     object= sxpar(header1, 'OBJECT')
     xrt_gti = mrdfits(filename_evt,2,/silent)
                                ;REMOVES SHORT GAPS IN GTIS
     gti_check, xrt_gti, xrt_gti_check
     swift_mkf = mrdfits(filename_mkf,1,/silent)
     swift_coordinates = mrdfits(filename_coord,2,/silent)
     swift_sao =  mrdfits(filename_sao,1,/silent)
     earth_angle = fltarr(n_elements(swift_sao.TIME))
     for k = 0l, n_elements(swift_sao.TIME) -1 do earth_angle(k) = Radiator_Earth_Angle(swift_sao(k).position, swift_sao(k).quaternion)*180./!PI
  endelse  
end
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;********************************
;START OF MAIN PROGRAM
;*******************************

pro brte_project_cla_earthangle, filename_evt, filename_mkf, filename_coord, filename_sao, g32

;COLORS:
DEVICE,DECOMPOSED=0.
tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1

;1;putting eventlist into structure in idl
print, 'Analysis of dataset: ',filename_evt
brte_read32, filename_evt, filename_mkf, filename_coord, filename_sao, g32, obsid, raobj, decobj, papnt, object, xrt_dataset, xrt_gti, swift_mkf, xrt_dataset_unfilt, stopper, swift_coordinates, swift_sao, earth_angle
if stopper eq 0 then begin ;if arrays are empty (ie, no bright earth) then DO NOTHING
    
;Define bright earth area coordinates:
   min_detx = 0
   max_detx = 80
   min_dety = 200
   max_dety = 400
   
;*** SELECT EVENTS FROM BRIGHT EARTH AREA
;xrt_dataset_croped_index = array indices that match the selected area
;xrt_dataset_crop = array of the selected area
;***
   cropping_data, xrt_dataset, min_detx, max_detx, min_dety, max_dety, xrt_dataset_crop, xrt_dataset_crop_2dhis, stopper
   
   if stopper eq 0 then begin ;if arrays are empty (ie, no bright earth) then DO NOTHING
;***DERIVE THRESHOLD FOR HOT PIXELS DEFINITION, REMOVE HOT PIXELS
;***NOTE: WE DON'T ACTUALLY USE THIS PART WHEN WE DO A GRADE32
;SELECTION, AS THAT EXCLUDES AUTOMATICALLY ALL THE HOT PIXELS
      threshold=threshold_deriver(xrt_dataset_crop_2dhis)
      index_match= exlude_hot_pixels(xrt_dataset_crop_2dhis, 100, xrt_dataset_crop)
;time_br =  xrt_dataset_crop[index_match].time
      time_br =  xrt_dataset_crop.time
                                ;***** PLOTS DATA FOR EACH ORBIT, SAVE INFO INTO FILES
      plotting_routine, obsid, xrt_gti, time_br, raobj, decobj, papnt, object, swift_mkf, xrt_dataset_unfilt, swift_coordinates, swift_sao, earth_angle
      
;***** WRITE FILE WITH CLEANED GTIs FOR EACH ORBITS (REMOVES GAPS)
      namedata= 'brightearth'+strtrim(string(obsid),2)+'.txt'
      openw,lun,namedata,/get_lun
      start_times = xrt_gti.start
      stop_times = xrt_gti.stop
      if n_elements(start_times) gt 1 then print_output = [transpose(start_times), transpose(stop_times)] else print_output = [start_times,stop_times]
      printf, lun, print_output
      Free_lun, lun
   endif
endif
end

