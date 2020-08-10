ima=readfits()
ima e' 200 x  200
cx cy   e' il centro che interessa e rr e' il raggio

mask=ima*0
for i=0,199
  for j=0,199
     if sqrt ((j-cx)^2+(j-cy)^2)  lt  rad then mask =1
endfor
endfor
iu=where(mask eq 1)
flux=total(ima[iu]
end
