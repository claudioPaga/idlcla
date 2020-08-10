function ploth_jak, x, bins, repeats

print, size(x)

y=dblarr(bins*repeats)
z=dblarr(bins*repeats)

h=histogram(x,binsize=(1D/double(bins)))
for i=0,(bins*repeats-1) do begin
    y[i]=h[i mod bins]
    z[i]=double(i)/bins
end 

plot,z,y/mean(1),psym=10,xtitle="phase",ytitle="counts"

oploterr,z,y/mean(1),sqrt(y)/mean(1)

openw,1,'ploth.dat'

for i=0,(bins*repeats-1) do begin
    printf,1,z[i],y[i],sqrt(y[i])
end 

close,1
return, y

end
