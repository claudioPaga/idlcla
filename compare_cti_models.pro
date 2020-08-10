pro compare_cti_models, qdpfile, ctifile1, ctifile2

  readcol, qdpfile, enqdp, binqdp, fluxqdp, errqdp, format='(d,d,d,d)'
  readcol, ctifile1, en1, bin1, flux1, err1, format='(d,d,d,d)'
  readcol, ctifile2, en2, bin2, flux2, err2, format='(d,d,d,d)'

  max_y = max([fluxqdp, flux1, flux2])

  DEVICE,DECOMPOSED=0.
  tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1

  plot, enqdp, fluxqdp, BACKGROUND=1,color=0, xtit = 'En (keV)', ytit = 'Cnts/s/keV', tit = 'Black = Original, Blue = Average CTI, Red = Pixel CTI', yr=[0,max_y] 
  oplot, en1, flux1, color=4
  oplot, en2, flux2, color=3

  r_split = strsplit(ctifile1, '.', /extract)
  entry_device = !d.name
  plotfilenameps = r_split[0]+'_compare_cti_models.ps'
  set_plot,'ps'
  device, filename = plotfilenameps, /color
  
  plot, enqdp, fluxqdp, BACKGROUND=1,color=0, xtit = 'En (keV)', ytit = 'Cnts/s/keV', tit = 'Black = Original, Blue = Average CTI, Red = Pixel CTI', yr=[0,max_y] 
  oplot, en1, flux1, color=4
  oplot, en2, flux2, color=3
  
  device,/close  
  set_plot,entry_device

end

