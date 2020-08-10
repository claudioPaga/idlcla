pro onlyplot_earth_angle_maxbrte, table_file, table_angle_elvmin, table_angle_peak, plotfile

readcol, table_file, obsid, orbit, orb_start, orb_end, roll_angle, sun_angle, min_elv, min_elv_with_sun, min_elv_with_sun_time, min_elv_brte_count, max_brte_count, time_of_max_brte, elv_at_max_brte, sun_when_max, long_at_start, lat_at_start, long_at_end, lat_at_end, format='(a,i,d,d,i,i,i,i,d,l,l,d,i,i,i,i,i,i)'
readcol, table_angle_elvmin, a1, a2, a3, a4, a5, a6, a7, a8, a9, earth_angle_array_minelv38, format = '(a, i, d, d, d, d, d, d, d, d)'
readcol, table_angle_elvmin, a1, a2, a3, a4, a5, a6, a7, a8, a9, earth_angle_array_max_brte, format = '(a, i, d, d, d, d, d, d, d, d)'

brte_cutoff = 50
bin_variable = 10.

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
DEVICE,DECOMPOSED=0.
tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1
!p.multi= [0,2,2,1,1]
nameps= 'earth_angle_hist.ps'
set_plot,'PS'
device,filename=nameps,    /color,/landscape ;, xs=30,ys=20

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
spawn,'gv earth_angle_hist.ps&'
;WRITE_PNG, 'earth_angle_brte_hist_plot_with_elv38.png', /transparent, tvrd(true =1)
end
