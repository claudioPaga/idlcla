pro lcomando,ra,dec
  evtfile=findfile('*.evt')
  tab=mrdfits(evtfile[0],1,hd1,/silent)
  trigstart=sxpar(hd1,'TSTART')
  trigstart=trigstart-1000.
  mkffile=findfile('*.mkf')
  satfile=findfile('*sat.fits')
  hkfile=findfile('*.hk')
  evtst=string(evtfile[0])
  mkfst=string(mkffile[0])
  satst=string(satfile[0])
  hkst=string(hkfile[0])
  rast=strtrim(string(ra),2)
  decst=strtrim(string(dec),2)
  frase='xrtgrblc clobber=yes evtfiles='+evtst+' mkffiles='+mkfst+' xhdfiles='+hkst+' attfiles='+satst+' srcra='+rast+' srcdec='+decst+ ' outstem=s1 pcsteprates=1:20,10:50,100:75,10000:150 wtsteprates=1:20,10:50,100:75,10000:150 splitorbits=yes bintype=1 cutlowbins=yes trigtime=2.484545457501600E+08 outdir=s1 minenergy=0.3 maxenergy=10.0 splitenergy=2.0'
  openw,lcom,'lcomando.txt',/get_lun
  printf,lcom,frase
  free_lun,lcom
  spawn,'source lcomando.txt'
end
