PRO bright_earth_angle , saafile, eclipsefile, earthfile, countfile, cnt_limit
;
;saafile=Text file with saa enter and saa exit times
;eclipsefile=Text file with eclipse in and out times
;earthfile=Text file with earth angles
;countfile=Text file with countrates
;modefiel=Text file with XRT mode (7=PC)
;eclipseflag:  'i'=Choose times while in eclipse   'o'=Choose times
;while out of eclipse
;
;This one gives the countrate when 


;Gets bad times for SAA 

readcol,saafile,saatimein,saatimeout,format='(a,a)' 
saain=strarr(n_elements(saatimein))
saaout=strarr(n_elements(saatimein))
for i=0,n_elements(saain)-1 do saain[i]=date2met(saatimein[i])  
for i=0,n_elements(saain)-1 do saaout[i]=date2met(saatimeout[i])   


;Gets bad times for Eclipse 

readcol,eclipsefile,ecltimein,ecltimeout,format='(a,a)' 
eclin=strarr(n_elements(ecltimein))
eclout=strarr(n_elements(ecltimein))
for i=0,n_elements(eclin)-1 do eclin[i]=date2met(ecltimein[i])  
for i=0,n_elements(eclin)-1 do eclout[i]=date2met(ecltimeout[i])   

;Reads elevation angles

readcol,earthfile,elvtime,elvang,format='(a,i)'     
elvtime_met=strarr(n_elements(elvtime))
for iii=0L,n_elements(elvtime)-1 do elvtime_met[iii]=date2met(elvtime[iii])

;Reads countrate 

readcol,countfile,cnttime,cnt,format='(a,f)'
cnttime_met=strarr(n_elements(cnttime))
for iiii=0L,n_elements(cnttime)-1 do cnttime_met[iiii]=date2met(cnttime[iiii])

elv_inter = INTERPOL(elvang,double(elvtime_met),double(cnttime_met),/SPLINE)
for cc=0,n_elements(elv_inter)-1 do begin
    if (elv_inter[cc] lt 0) then elv_inter[cc]=elv_inter[cc-1]
endfor

print,elv_inter
elv_limit=fltarr(n_elements(elv_inter))

count_index=0
while count_index lt n_elements(cnt)-1 do begin
    
    count_index=count_index-1
    repeat count_index=count_index+1 until ((cnt[count_index] gt cnt_limit) or (count_index eq n_elements(cnt)-1))
    time_cnt=cnttime_met[count_index]
   ; print,'Indice ',count_index
   ; print,'Conteggi: ',cnt[count_index]
   ; stop

;Check if time is ok (no CCD Temp, SAA)
    
    flag=0
    
    for check=0,n_elements(saain)-1 do begin ;check for SAA
        if ((time_cnt gt saain[check]) and (time_cnt lt saaout[check])) then flag=1
    endfor 
    
    
    if flag eq 0 then begin     ;If not SAA then check for eclipse
        for check_ecl=0,n_elements(eclin)-1 do begin ;check for Eclipse
            if ((time_cnt gt eclin[check_ecl]) and (time_cnt lt eclout[check_ecl])) then flag=1 ;If in eclipse ===> flag=1
        endfor 
    endif
        
    if flag eq 0 then begin
        elv_limit[count_index]=elv_inter[count_index]        
    endif

    print,'Elv angle for crossing the countlimit: ',elv_limit[count_index]
   ; stop
    while (cnt[count_index] gt cnt_limit) do count_index=count_index+1 
   ; print,'count index ',count_index
   ; stop
endwhile
;print,elv_limit
index=where(elv_limit, num)
elv_limit_corr=fltarr(num)
elv_limit_corr=elv_limit[index]

print,elv_limit_corr

printname='elv_limit.txt'

openw,nomelogico,printname,/get_lun
printf,nomelogico,elv_limit_corr
free_lun,nomelogico

angltito=strmid(string(cnt_limit),4,4)

xtito='Earth Elevation angle'
tito='Elvation angle for countrate > '+angltito+' cnts/s'

set_plot,'X'

plothisto,elv_limit_corr,100,160,5,xtito,tito

psname='elv_limit_'+angltito+'.ps'

set_plot,'PS'
device,filename=psname;,xs=35,ys=30,/col

plothisto,elv_limit_corr,100,160,5,xtito,tito

device,/close


end
