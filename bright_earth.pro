PRO bright_earth , saafile, eclipsefile, earthfile, countfile, modefile, eclipseflag
;
;saafile=Text file with saa enter and saa exit times
;eclipsefile=Text file with eclipse in and out times
;earthfile=Text file with earth angles
;countfile=Text file with countrates
;modefiel=Text file with XRT mode (7=PC)
;eclipseflag:  'i'=Choose times while in eclipse   'o'=Choose times
;while out of eclipse



;Gets bad times for SAA 

readcol,saafile,saatimein,saatimeout,format='(a,a)' 
saain=strarr(n_elements(saatimein))
saaout=strarr(n_elements(saatimein))
for i=0,n_elements(saain)-1 do saain[i]=date2met(saatimein[i])  
for i=0,n_elements(saain)-1 do saaout[i]=date2met(saatimeout[i])   

;Gets good times for PC mode

readcol,modefile,modetime,modenum,format='(a,a)'
k=0
for j=0,n_elements(modetime)-1 do if (modenum[j] eq '7') then k=k+1
modetimein=strarr(k)
modetimeout=strarr(k)
k=0
for j=0,n_elements(modetime)-1 do begin
    if (modenum[j] eq '7') then begin
        modetimein[k]=modetime[j] 
        modetimeout[k]=modetime[j+1]  
        k=k+1  
    endif
endfor

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

;elv_range=[97,101,103,105,107,109,111,113,115,117,119]

elv_range=[97,105,115,125,135,145,155]


cntang_range=fltarr(n_elements(elv_range)-1,201)

;I'll go trough each elv_angles ranges with a for loop

for elv_co=0,n_elements(elv_range)-2 do begin
   
    elv_min=elv_range[elv_co]
    elv_max=elv_range[elv_co+1]
    print,'Elevation angle range: ',elv_min,elv_max

;Second for loop to go trough all the elv angles
    for elv_index=0L,n_elements(elvang)-1 do begin
        if ((elvang[elv_index] ge elv_min) and (elvang[elv_index] lt elv_max)) then begin
            
            elv_time=elvtime_met[elv_index]

;Check if time is ok (no CCD Temp, SAA, mode constraint)
            
            flag=0

            for check=0,n_elements(saain)-1 do begin ;check for SAA
                if ((elv_time gt saain[check]) and (elv_time lt saaout[check])) then flag=1
            endfor 
            
            
            if flag eq 0 then begin ;If not SAA then check for PC modes
                
                for cont=0,n_elements(modetimein)-1 do begin
                    if ((elvtime[elv_index] gt modetimein[cont]) and (elvtime[elv_index] lt modetimeout[cont])) then flag=2
                endfor
            endif
            
            if flag eq 2 then begin ;If PC mode then check for Eclipse
                
                for check_ecl=0,n_elements(eclin)-1 do begin ;check for Eclipse
                    if ((elv_time gt eclin[check_ecl]) and (elv_time lt eclout[check_ecl])) then flag=3     ;If in eclipse ===> flag=3
                endfor 
            endif
            
            
            
            if (eclipseflag eq 'i' and flag eq 3) then flag=4
            if (eclipseflag eq 'o' and flag eq 2) then flag=4
                        
            if flag eq 4 then begin
                
                time_n=0L ;Look what the coutrate is for that elv angle time
                repeat begin
                    timects=cnttime_met[time_n]
                    time_n=time_n+1
                endrep until (timects gt elv_time)
                if (cnt[time_n-1] le 200) then begin
                    cnt_round=round(cnt[time_n-1])
                    cntang_range[elv_co,cnt_round]= cntang_range[elv_co,cnt_round]+1 
                endif else  cntang_range[elv_co,200]=cntang_range[elv_co,200]+1
            endif
        endif
    endfor
endfor



cntang_range_tr=transpose(cntang_range)
print,cntang_range_tr

if eclipseflag eq 'i' then printname='elv_table_inecl.txt'
if eclipseflag eq 'o' then printname='elv_table_outecl.txt'


openw,nomelogico,printname,/get_lun
printf,nomelogico,cntang_range_tr
free_lun,nomelogico


total_cnts=intarr(n_elements(elv_range)-1)
;cts_range=[0,5,10,20,50,99,101]

cts_range=[0,10,20,50,200]
cts_range_tot=intarr(n_elements(elv_range)-1,n_elements(cts_range)-1)
cts_range_per=fltarr(n_elements(elv_range)-1,n_elements(cts_range)-1)


for elv_cc=0,n_elements(elv_range)-2 do begin
    total_cnts[elv_cc]=total(cntang_range[elv_cc,*])
    for aa=0,n_elements(cts_range)-2 do begin
        for dd=(cts_range[aa]),(cts_range[aa+1]-1) do begin
            cts_range_tot[elv_cc,aa]=cts_range_tot[elv_cc,aa]+cntang_range[elv_cc,dd]
        endfor
        cts_range_per[elv_cc,aa]=float(cts_range_tot[elv_cc,aa])/float(total_cnts[elv_cc])
    endfor 
endfor


for ii=0,n_elements(elv_range)-2 do begin
    elv_min=elv_range[ii]
    elv_max=elv_range[ii+1]
    print,'Elevation angle range: [',elv_min,elv_max,']'
    print,cts_range_per[ii,*]
endfor

cts_range_per_tr=transpose(cts_range_per)
print,cts_range_per
openw,logico,'elv_histo_allmodes_ecl_very.txt',/get_lun
printf,logico,cts_range_per_tr
free_lun,logico

;xplot_range=[99,102,104,106,108,110,112,114,116,118]

xplot_range=[100,110,120,135,140,150]

if eclipseflag eq 'i' then psname='bright_earth_inecl.ps'
if eclipseflag eq 'o' then psname='bright_earth_outecl.ps'

set_plot,'PS'
device,filename=psname;,xs=35,ys=30,/col

!p.multi=[0,3,2,0,0]
for plot_c=0,n_elements(cts_range)-2 do begin
    r1=strmid(string(cts_range[plot_c]),5,4)
    r2=strmid(string(cts_range[plot_c+1]),5,4)
    if eclipseflag eq 'i' then titolo='Counts ['+r1+r2+']   IN ECLIPSE'
    if eclipseflag eq 'o' then titolo='Counts ['+r1+r2+']   OUT OF ECLIPSE'
    plot,xplot_range,cts_range_per[*,plot_c],psym=4,yr=[0,1],tit=titolo,xtit='Elev angle',ytit='Cnts %',background=16777215,color=0
endfor
!p.multi=0
device,/close


;x_plot_cts=[2.5,7.5,15,25,75,100]

x_plot_cts=[5,15,25,75,200]


set_plot,'PS'

if eclipseflag eq 'i' then psname='bright_earth_elv_inecl.ps'
if eclipseflag eq 'o' then psname='bright_earth_elv_outecl.ps'


device,filename=psname;,xs=35,ys=30,/col
!p.multi=[0,4,3,0,0]
for plot_b=0,n_elements(xplot_range)-1 do begin
    r1=strmid(string(elv_range[plot_b]),5,4)
    r2=strmid(string(elv_range[plot_b+1]),5,4)
    if eclipseflag eq 'i' then titolo='Elv ang ['+r1+r2+']   IN ECLIPSE'
    if eclipseflag eq 'o' then titolo='Elv ang ['+r1+r2+']   OUT OF ECLIPSE'
    plot,x_plot_cts,cts_range_per[plot_b,*],psym=4,yr=[0,1],tit=titolo,xtit='Cnts range',ytit='Cnts %',background=16777215,color=0
endfor
!p.multi=0


device,/close




end
