pro plot_grbcla,fname,outname,x1,x2,y1,y2,f1,f2
!p.multi=0
setup_colors
openr,lun,fname,/get_lun
set_plot,'ps'
device,filename=outname,/color
line=parse_str(lun)
line=parse_str(lun)
line=parse_str(lun)
line=float(line)
if line(10) ge 3.0 then begin
    if line(13) eq 0 then plot,[line(1),line(2)],[line(3),line(3)],/xlog,/ylog,xrange=[x1,x2],yrange=[y1,y2],xtitle='Seconds since BAT trigger',ytitle='Counts/s',tit=outname,charsize=1.5,xstyle=1,ystyle=1
    if line(13) eq 1 then plot,[line(1),line(2)],[line(3),line(3)],/xlog,/ylog,xrange=[x1,x2],yrange=[y1,y2],xtitle='Seconds since BAT trigger',ytitle='Counts/s',tit=outname,charsize=1.5,xstyle=1,ystyle=1
    if line(13) eq 0 then oplot,[line(1),line(2)],[line(3),line(3)],color=100
    if line(13) eq 1 then oplot,[line(1),line(2)],[line(3),line(3)],color=200
    if line(13) eq 0 then oplot,[line(0),line(0)],[line(3)-line(4),line(3)+line(4)],color=100
    if line(13) eq 1 then oplot,[line(0),line(0)],[line(3)-line(4),line(3)+line(4)],color=200
endif
if line(10) lt 3.0 then begin
    if line(13) eq 0 then plot,[line(1),line(2)],[line(3),line(3)],/xlog,/ylog,xrange=[x1,x2],yrange=[y1,y2],xtitle='Seconds since BAT trigger',ytitle='Counts/s',tit=outname,charsize=1.5,xstyle=1,ystyle=1
    if line(13) eq 1 then plot,[line(1),line(2)],[line(3),line(3)],/xlog,/ylog,xrange=[x1,x2],yrange=[y1,y2],xtitle='Seconds since BAT trigger',ytitle='Counts/s',ctit=outname,harsize=1.5,xstyle=1,ystyle=1
    if line(13) eq 0 then oplot,[line(1),line(2)],[line(3)+line(4),line(3)+line(4)],color=100
    if line(13) eq 1 then oplot,[line(1),line(2)],[line(3)+line(4),line(3)+line(4)],color=200
    if line(13) eq 0 then arrow,line(0),line(3)+line(4),line(0),line(3)-line(4),/data,color=100
    if line(13) eq 1 then arrow,line(0),line(3)+line(4),line(0),line(3)-line(4),/data,color=200
endif

while not eof(lun) do begin 
    line=parse_str(lun)
    line=float(line)
    if line(10) ge 3.0 then begin
        if line(13) eq 0 then oplot,[line(1),line(2)],[line(3),line(3)],color=100
        if line(13) eq 0 then oplot,[line(0),line(0)],[line(3)-line(4),line(3)+line(4)],color=100
        if line(13) eq 1 then oplot,[line(1),line(2)],[line(3),line(3)],color=200
        if line(13) eq 1 then oplot,[line(0),line(0)],[line(3)-line(4),line(3)+line(4)],color=200
    endif
    if line(10) lt 3.0 then begin
        if line(13) eq 0 then oplot,[line(1),line(2)],[line(3)+line(4),line(3)+line(4)],color=100
        if line(13) eq 1 then oplot,[line(1),line(2)],[line(3)+line(4),line(3)+line(4)],color=200
        if line(13) eq 0 then arrow,line(0),line(3)+line(4),line(0),line(3)-line(4),/data,color=100
        if line(13) eq 1 then arrow,line(0),line(3)+line(4),line(0),line(3)-line(4),/data,color=200
    endif
endwhile

risposta=''
read,risposta,prompt='Overplot? [Y/N]'
if risposta eq 'Y' or risposta eq 'y' then begin

    xp=indgen(1000)/100.0+4.5
    yp=f1+f2*xp
    xpt=10^xp
    ypt=10^yp
    oplot,xpt,ypt
endif
    free_lun,lun


device,/close
end
