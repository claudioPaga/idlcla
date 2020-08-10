pro rayleighs_pc, ev, period, step, number, n

; Do the Z^2_n test of Bucherri et al. (1983)

pc=2.507

zsquare = dblarr(number)

xaxis = dblarr(number)

p = period - double(number/2)*step

s=long(size(ev))
points=s[1]

seed = 1234

;for i=0L,points-1 do begin
;    ev[i].barytime = ev[i].barytime+(pc*randomu(seed))
;end 

pi = 2D*acos(0D)

print,'points =', points

for i=0L, number-1 do begin
    phase = fold(ev.barytime,p,ev[0].time)
;    plot,lc
    for k=1L, n do begin
        coses = 0
        sines = 0
        for j=0L, points-1 do begin
            coses = coses+cos(double(k)*phase[j]*2D*pi)
            sines = sines+sin(double(k)*phase[j]*2D*pi)
        end
;        print, 'k =',k
        zsquare[i] = zsquare[i] + coses^2D + sines^2D
    end
    zsquare[i] = zsquare[i] *(2.0D/points)
    xaxis[i]=p-period
    p = p + step
end


; Set up postscript file plot
;set_plot,'ps'
;device, filename='rayleigh.ps', /landscape
xper = xaxis + period
; plot the periodogram


plot, xper, zsquare,yrange=[0,max(zsquare)],xrange=[min(xper),max(xper)],xstyle=1, $
  xtitle = "Period (s)",ytitle="Z!S!E2!N!R!I1",xtickformat='(F9.7)',xticks=6
oplot,[min(xper),max(xper)],[15.9361,15.9361],linestyle=1
xyouts,min(xper)+20*step,15.9361+0.1,"90% Confidence"
oplot,[min(xper),max(xper)],[21.286,21.286],linestyle=1
xyouts,min(xper)+20*step,21.286+0.1,"99% Confidence"
oplot,[period,period],[0,30],linestyle=2
;close the file
;device,/close
set_plot,'x'


peak = xaxis(where(zsquare eq max(zsquare)))

print,"delta P = ",peak,format='(A,E20.9)'
print,"Period  = ",peak+period,format='(A,D20.11)'
print,"Step    = ",step,format='(A,D20.11)' 
print,"Height  = ",max(zsquare),format='(A,D20.11)' 

openw,1,'rayleighs.out'
for i = 0, number-1 do begin
    printf,1,xaxis[i],zsquare[i]
end

close,1
end

    
