pro silicon_ccd_histo,filename
DEVICE,DECOMPOSED=0.
tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1


;PURPOSE: Create histogram of the number of events in the 
tab=mrdfits(filename,1,hd1)

indexsi = where (tab.PI gt 170 and tab.PI lt 200 and tab.DETX ge 200 and tab.DETX lt 400 and tab.DETY ge 200 and tab.DETY lt 400, n_indexsi)
print,'Total number of Si events = ',n_indexsi

x_silicon=tab(indexsi).DETX
y_silicon=tab(indexsi).DETY
n_si_events=intarr(200,200)

for index_x=200,399 do begin
    
    for index_y=200,399 do begin
        xyindex = where (x_silicon eq index_x and y_silicon eq index_y, xynumber)
        n_si_events(index_x-200,index_y-200) = xynumber
    endfor
endfor

;here i "take out" the hot columns that would give an incorrect high number of 

n_si_events(290-200,*)=-100
n_si_events(291-200,*)=-100
n_si_events(292-200,*)=-100
n_si_events(293-200,*)=-100
n_si_events(294-200,*)=-100
;n_si_events(295-200,*)=-100
n_si_events(319-200,*)=-100
n_si_events(320-200,*)=-100
n_si_events(321-200,*)=-100
index_gt10 = where (n_si_events gt 10, n_index_gt10)
print,'Number of pixels in the central CCD with > 10 Si events:',n_index_gt10
index_eq0 = where (n_si_events eq 0, n_index_eq0)
print,'Number of pixels in the central CCD with 0 Si events:',n_index_eq0

maxx = fix(max (n_si_events) + max (n_si_events)*0.1)

plothist, n_si_events, xhist, yhist, color=0, background=1, /fill, fcolor=2, xr=[-1,maxx], xsty=1, xtit='Number of events',tit='Distribution of Silicon events per pixel',ytit='Number of pixels'

positiv_index =  where (xhist ge 0)
tot_eve = total (yhist(positiv_index))
aaa=[transpose(xhist(positiv_index)-0.5),transpose(yhist(positiv_index)*1./tot_eve*1.)]

print, aaa
WRITE_PNG, 'histo_si.png', TVRD(TRUE=1), /transparent
end
