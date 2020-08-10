pro rayleighs2, time, period, step, number, n

; Do the Z^2_n test of Bucherri et al. (1983)

zsquare = dblarr(number)

xaxis = dblarr(number)

p = period - double(number/2)*step

s=size(time)
points=s[1]

pi = 2D*acos(0D)

print,'points =', points

for i=0L, number-1 do begin
    phase = fold(time,p,time[0])
    for k=1, n do begin
        coses = 0
        sines = 0
        for j=0L, points-1 do begin
            coses = coses+cos(double(k)*phase[j]*2D*pi)
            sines = sines+sin(double(k)*phase[j]*2D*pi)
        end
        zsquare[i] = zsquare[i] + coses^2D + sines^2D
    end
    zsquare[i] = zsquare[i] *(2.0D/points)
    xaxis[i]=p-period
    p = p + step
end

; Find best fit period
peak = xaxis(where(zsquare eq max(zsquare)))

P = peak + period

; Plot periodogram and folded light curve
window,0
plot, xaxis, zsquare
window,1
print,P[0]
x=ploth_jak(fold(time,P[0],time[0]),16,3)


; Print out best fit period, delta P and error
print,"delta P = ",peak,format='(A,E20.9)'
print,"Period  = ",peak+period,format='(A,D20.11)'
print,"Error   = ",step,format='(A,D20.11)'

end

    
