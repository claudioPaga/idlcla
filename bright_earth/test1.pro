pro test1

;this program will find brocken detector cells of a bright earth image and remove their reading from an event list

;COLORS:

;(4 colors)

DEVICE,DECOMPOSED=0.
tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1

;plot,x,y,BACKGROUND=1,color=0 (white BG, black points)

;color=0 black
;color=2 green
;color=3 red
;color=4 blu



;1;putting eventlist into structure in idl
xrt_dataset = mrdfits('sw00528195000xpcw3po_uf.evt',1, header1)

obsid= sxpar(header1, 'OBS_ID')
raobj= sxpar(header1, 'RA_OBJ')
decobj= sxpar(header1, 'DEC_OBJ')
papnt= sxpar(header1, 'PA_PNT')
object= sxpar(header1, 'OBJECT')

help, /str, xrt_dataset


xrt_gti = mrdfits('sw00528195000xpcw3po_uf.evt',2)


swift_attitude = mrdfits('sw00528195000s.mkf',1)


;making brite earth contaminated region array, xrt_dataset only contains
;array of indices

;Define bright earth area coordinates:
;
min_detx = 0
max_detx = 80
min_dety = 200
max_dety = 400

;2;Filter on bright earth area
;xrt_dataset_croped_index = array indices that match the selected area
;xrt_dataset_crop = array of the selected area

xrt_dataset_croped_index = where(xrt_dataset.detx gt min_detx and xrt_dataset.detx lt max_detx and xrt_dataset.dety gt min_dety and xrt_dataset.dety lt max_dety, total_matches)
xrt_dataset_crop = xrt_dataset(xrt_dataset_croped_index)

;plotting the histogram
plothist, xrt_dataset_crop.time

;Bright earth vs hot pixels selection

;3;making 2d hist
;We need the _sub array because HIST_2D returns an array that is
;max_detx x max_dety long
xrt_dataset_almostcrop_2dhis= HIST_2D( xrt_dataset_crop.detx, xrt_dataset_crop.dety)
xrt_dataset_crop_2dhis = xrt_dataset_almostcrop_2dhis[min_detx:max_detx-1, min_dety:max_dety-1]

;4;findig total 0's in sub
xrt_dataset_crop_his_no0_index = where(xrt_dataset_crop_2dhis ne 0)

xrt_dataset_crop_his_no0 = xrt_dataset_crop_2dhis(xrt_dataset_crop_his_no0_index)


median_be = MEDIAN (xrt_dataset_crop_his_no0)

variance_be = VARIANCE (xrt_dataset_crop_his_no0)

y_xrt_dataset_crop_his_no0_his= histogram(xrt_dataset_crop_his_no0, bin=1, location= x_numberofphotonhits)



;4.5; method 1 for finding threshhold
;assuming pixel distribution is Gaussian,
;sort the numbers

xrt_dataset_crop_his_no0_index_sorted= sort(xrt_dataset_crop_his_no0)
xrt_dataset_crop_his_no0_sorted=xrt_dataset_crop_his_no0[xrt_dataset_crop_his_no0_index_sorted]


xrt_dataset_crop_his_no0_sorted_66index = float(n_elements(xrt_dataset_crop_his_no0_sorted)*(2./3.))
xrt_dataset_crop_his_no0_sorted_66 = xrt_dataset_crop_his_no0_sorted[fix(xrt_dataset_crop_his_no0_sorted_66index)]
threshold_3sig= 3*xrt_dataset_crop_his_no0_sorted_66


;4.6
brte_index = where(xrt_dataset_crop_2dhis le threshold_3sig)

brte_index_2d =  ARRAY_INDICES(xrt_dataset_crop_2dhis, brte_index)

brte_index_2d[1,*] = brte_index_2d[1,*]+200



length= fix(n_elements(brte_index_2d)/2.)


time_br = [0]


xrt_dataset_string =  strtrim(string(xrt_dataset_crop.detx), 2)+strtrim(string(xrt_dataset_crop.dety), 2)

brte_index_string =  strtrim(string(brte_index_2d[0,*]), 2 )+strtrim(string(brte_index_2d[1,*]), 2)

brte_index_string1d =  strarr(n_elements(brte_index_string))

for i =  0,  n_elements(brte_index_string)-1 do brte_index_string1d[i] =   brte_index_string[i]

match2, xrt_dataset_string,  brte_index_string1d,  sub1,  sub2

indexmatch =  where(sub1 ge 0)

time_br =  xrt_dataset_crop[indexmatch].time

stop

for i= 1 ,length-1 do begin




    brte_index=where(xrt_dataset_crop.detx eq brte_index_2d[0,i] and xrt_dataset_crop.dety eq brte_index_2d[1,i], n_brte_index)
    if n_brte_index ge 1 then time_br =[time_br,xrt_dataset_crop[brte_index].time]


;    if n_brte_index eq 1 and brte_index ne -1 then time_br =[time_br,xrt_dataset_crop[brte_index].time]
  
endfor
time_br = time_br[1:n_elements(time_br)-1]
plothist, time_br



;plotting routine

nomeps= 'plot1.ps'
set_plot,'PS'
device,filename=nomeps,    /color, xs=30,ys=20











!p.multi= [0,4,2,0,1]
for j= 0,n_elements(xrt_gti.start)-1 do begin
    histogram_time_br= histogram( time_br, bin=20)
    plothist, time_br-xrt_gti[j].start, xr=[-10,xrt_gti[j].stop-xrt_gti[j].start],/ylog, yr=[0.1, max(histogram_time_br)], bin=20,xtit='time', ytit='counts', tit= 'bright earth count rate', background= 1, color=0
    
    screen_text1 = 'obsid = '+strtrim(string(obsid),2)
    XYOUTS, 0.65, 0.9, screen_text1, /norm
    screen_text = '(RA, DEC) = '+strtrim(string(raobj),2)+','+strtrim(string(decobj),2)
    XYOUTS, 0.65, 0.85, screen_text, /norm

    screen_text2 = 'roll angle = '+strtrim(string(papnt),2)
    XYOUTS, 0.65, 0.8, screen_text2, /norm

    screen_text3 = 'object name = '+strtrim(string(object),2)
    XYOUTS, 0.65, 0.75, screen_text3, /norm

orbit_index = where(swift_attitude.time ge xrt_gti[j].start and swift_attitude.time le xrt_gti[j].stop)
orbit= swift_attitude[orbit_index]
orbit_rollangle = median(orbit.roll)

screen_orbit= 'orbit_rollangle = '+strtrim(string(orbit_rollangle),2)

plot, swift_attitude.time-xrt_gti[j].start, swift_attitude.elv,  xr=[-10,xrt_gti[j].stop-xrt_gti[j].start], yr=[0,100], color=0, title=screen_orbit


swift_attitude_sun_index = where(swift_attitude.sunshine eq 1)
oplot, swift_attitude[swift_attitude_sun_index].time-xrt_gti[j].start, swift_attitude[swift_attitude_sun_index].elv, color=3, psym=3




endfor

device,/close
spawn,'gv '+nomeps+'&'

set_plot,'x'


end
