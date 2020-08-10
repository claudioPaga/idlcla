pro plot_roll_maxbrte, table_file

DEVICE,DECOMPOSED=0.
tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1


readcol, table_file, obsid, orbit, orb_start, orb_end, roll_angle, sun_angle, min_elv, min_elv_with_sun, min_elv_with_sun_time, min_elv_brte_count, max_brte_count, time_of_max_brte, elv_at_max_brte, sun_when_max, long_at_start, lat_at_start, long_at_end, lat_at_end, format='(a,i,d,d,i,i,i,i,d,l,l,d,i,i,i,i,i,i)'

;roll_angle = table_file.roll_angle
;max_brte_count = table_file.max_brte_count

brte_cutoff = 50

brte_orbs =  where(max_brte_count ge brte_cutoff)
brte_roll_angle = roll_angle[brte_orbs]
brte_roll_angle = [0,brte_roll_angle]

equlipse_at_min_elv = where(min_elv[brte_orbs] ne min_elv_with_sun[brte_orbs])
equlipse_roll_angle = roll_angle[equlipse_at_min_elv]

bin_variable = 10.

elv38_index = where(min_elv_with_sun le 38)
elv38_roll_angle = roll_angle[elv38_index]
max_brte_count_elv38 =  min_elv_brte_count[elv38_index]
brte_orbs_elv38 =  where(max_brte_count_elv38 ge brte_cutoff)
brte_elv38_roll_angle = elv38_roll_angle[brte_orbs_elv38]
brte_elv38_roll_angle = [0,brte_elv38_roll_angle]

max_brte_count_elv38_fullorb =  max_brte_count[elv38_index]
brte_orbs_elv38_fullorb =  where(max_brte_count_elv38_fullorb ge brte_cutoff)
brte_elv38_roll_angle_fullorb = elv38_roll_angle[brte_orbs_elv38_fullorb]

!p.multi= [0,1,2]

plothist, roll_angle, xroll, yroll, background = 1, color = 0, ytit='n orbits', xtit= 'roll angle', tit = 'bright earth roll angle distribution', bin = bin_variable, xr = [0,360], xstyle = 1, yr = [0,4000]
plothist, brte_roll_angle, xroll_brte, yroll_brte, /overplot, color = 3, /fill, fcolor = 3, bin = bin_variable

max_yroll = max(yroll)

roll_histo = histogram(roll_angle, bin = bin_variable, min = 0, locations = xrollhisto)
roll_brte_histo = histogram(brte_roll_angle, bin = bin_variable, min = 0, locations = xrollhisto_brte)

;plot, xroll, yroll_brte*1./ yroll, background = 1, color = 0,psym=2
;stop

plot, xrollhisto, roll_brte_histo*1./roll_histo, background = 1, color = 0, psym=2, tit = 'fraction of orbits affected by bright earth', xtit = 'roll angle', ytit = 'fraction affected by bright earth',xstyle = 1

screen_text1 = 'total number of orbits'
XYOUTS, 0.55, 0.9, screen_text1, /norm, color = 0, charsize = 1.2
screen_text2 = 'orbits affected by bright earth'
XYOUTS, 0.55, 0.87, screen_text2, /norm, color = 3, charsize = 1.2


WRITE_PNG, 'roll_angle_brte_hist_plot.png', /transparent, tvrd(true =1)

!p.multi= [0,3,2,1,1]
;**************************


plothist, roll_angle, xroll, yroll, background = 1, color = 0, ytit='n orbits', xtit= 'roll angle', bin = bin_variable, xr = [0,360], xstyle = 1, yr = [0,max_yroll], charsize = 1.3
plothist, brte_roll_angle, xroll_brte, yroll_brte, /overplot, color = 3, /fill, fcolor = 3, bin = bin_variable
;plothist, equlipse_roll_angle, xroll_brte, yroll_brte, /overplot, color = 0, /fill, fcolor = 2, bin = bin_variable
roll_histo = histogram(roll_angle, bin = bin_variable, min = 0, locations = xrollhisto)
roll_brte_histo = histogram(brte_roll_angle, bin = bin_variable, min = 0, locations = xrollhisto_brte)


;plot, xroll, yroll_brte*1./ yroll, background = 1, color = 0,psym=2
;stop

plot, xrollhisto, roll_brte_histo*1./roll_histo, background = 1, color = 0, psym=2, xtit = 'roll angle', ytit = 'fraction affected by bright earth',xstyle = 1, charsize = 1.3




plothist, elv38_roll_angle, xroll38, yroll38, background = 1, color = 0, ytit='n orbits', xtit= 'roll angle', tit = 'bright earth roll angle distribution', bin = bin_variable, xr = [0,360],xstyle = 1, yr = [0,max_yroll], charsize = 1.3
plothist, brte_elv38_roll_angle, xroll_brte38, yroll_brte38, /overplot, color = 3, /fill, fcolor = 3, bin = bin_variable

roll_38_histo = histogram(elv38_roll_angle, bin = bin_variable, min = 0, locations = xrollhisto_38)
roll_brte_38_histo = histogram(brte_elv38_roll_angle, bin = bin_variable, min = 0, locations = xrollhisto_38_brte)


;plot, xroll, yroll_brte*1./ yroll, background = 1, color = 0,psym=2
;stop

plot, xrollhisto_38, roll_brte_38_histo*1./roll_38_histo, background = 1, color = 0, psym=2, tit = 'fraction of orbits affected by bright earth', xtit = 'roll angle', ytit = 'fraction affected by bright earth',xstyle = 1, charsize = 1.3, xr = [0,360], yr = [0,1]

brte_elv38_roll_angle_fullorb = [0, brte_elv38_roll_angle_fullorb]

plothist, elv38_roll_angle, xroll38_fullorb, yroll38_fullorb, background = 1, color = 0, ytit='n orbits', xtit= 'roll angle', bin = bin_variable, xr = [0,360],xstyle = 1, yr = [0,max_yroll], charsize = 1.3
plothist, brte_elv38_roll_angle_fullorb, xroll_brte38_fullorb, yroll_brte38_fullorb, /overplot, color = 3, /fill, fcolor = 3, bin = bin_variable

roll_38_histo_fullorb = histogram(elv38_roll_angle, bin = bin_variable, min = 0, locations = xrollhisto_38_fullorb)
roll_brte_38_histo_fullorb = histogram(brte_elv38_roll_angle_fullorb, bin = bin_variable, min = 0, locations = xrollhisto_38_brte_fullorb)


;plot, xsun, ysun_brte*1./ ysun, background = 1, color = 0,psym=2
;stop

plot, xrollhisto_38_fullorb, roll_brte_38_histo_fullorb*1./roll_38_histo_fullorb, background = 1, color = 0, psym=2, xtit = 'roll angle', ytit = 'fraction affected by bright earth',xstyle = 1, xr = [0,360], yr = [0,1], charsize = 1.3


screen_text1 = 'all orbits'
XYOUTS, 0.07, 0.95, screen_text1, /norm, color = 0, charsize = 0.8
screen_text2 = 'all orbits with bright earth(BE)'
XYOUTS, 0.07, 0.93, screen_text2, /norm, color = 3, charsize = 0.8
;screen_text3 = 'orbits in eclipse'
;XYOUTS, 0.4, 0.835, screen_text3, /norm, color = 2, charsize = 0.8
screen_text1 = 'sunlit orbits (elv<38)'
XYOUTS, 0.41, 0.95, screen_text1, /norm, color = 0, charsize = 0.8
screen_text2 = 'sunlit orbits with BE @ elv<38'
XYOUTS, 0.41, 0.93, screen_text2, /norm, color = 3, charsize = 0.8
screen_text1 = 'sunlit orbits (elv<38)'
XYOUTS, 0.74, 0.95, screen_text1, /norm, color = 0, charsize = 0.8
screen_text2 = 'sunlit orbits (elv<38) with BE'
XYOUTS, 0.74, 0.93, screen_text2, /norm, color = 3, charsize = 0.8

WRITE_PNG, 'roll_angle_brte_hist_plot_with_elv38.png', /transparent, tvrd(true =1)

;%%%%%%%%%%%%%%%%%%%%%%%%%%new%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!p.multi= [0,4,2,0,0]
nameps= 'elv_vs_roll_sets.ps'
set_plot,'PS'
device,filename=nameps,    /color,/landscape

for j = 0, (360/bin_variable) - 1 do begin
set_roll_index = where(roll_angle ge j*bin_variable and roll_angle le (j*bin_variable)+bin_variable and max_brte_count ge brte_cutoff, total)
if total gt 1 then begin 
set_roll_elv  = elv_at_max_brte[set_roll_index]+ randomn(seed, n_elements(elv_at_max_brte[set_roll_index]))

graph_title = 'brte elv at roll ='+strtrim(string(fix(j*bin_variable)),2)+'-'+strtrim(string(fix((j+1)*bin_variable)),2)
plothist, set_roll_elv, bin = 1, xr = [25,60], xstyle = 1, background = 1, color = 0, ytit='n orbits', xtit= 'elv angle', tit = graph_title, charsize = 1.3
endif
endfor

device,/close



!p.multi= 0
nameps= 'elv_vs_roll.ps'
set_plot,'PS'
device,filename=nameps,    /color,/landscape
plot, roll_angle[brte_orbs]+ randomn(seed, n_elements(elv_at_max_brte[brte_orbs])), elv_at_max_brte[brte_orbs]+ randomn(seed, n_elements(elv_at_max_brte[brte_orbs])), xr = [0,360], psym = 3, background = 1, color = 0, tit = 'elv at bright earth peak', ytit = 'elv', xtit = 'roll'

device,/close

end
