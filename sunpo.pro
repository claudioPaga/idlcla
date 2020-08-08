pro sunpo

anno=2006
read,anno,prompt='Year:(2006)'
mese=0
giorno=0
read,mese,prompt='Month (mm):'
read,giorno,prompt='Day (dd):'
  jdcnv, anno, mese, giorno, 0 ,jd      ;Find Julian date jd = 2445090.5   
  sunpos, jd, ra, dec
  print,adstring(ra,dec,2)
  rad=TEN(ra)
  decd=TEN(dec)
  print,rad,decd
end

