pro phil2spec,inputlc
  
                                ;takes lc from phil's leicester website
                                ;modifies it so it can be read by lorella's format.pl script
  
  ;INPUT: phil_lcname, ex: 'lc.txt'
  
  readcol,inputlc,t,tp,tm,cr,crerr,format='(d,d,d,d,d)'
  t1=t+tm
  t2=tp-tm
  openw,nologico,'lc_lorella.txt',/get_lun
  outvec=[transpose(t1),transpose(t2),transpose(cr),transpose(crerr)]
  printf,nologico,outvec
  free_lun,nologico
  
  spawn,'~/format2.pl -i lc_lorella.txt -o lcout.dat -r lc_xspec'
  tab=mrdfits('lc_xspec.pha',1,hd1)
  print,'xspec'
  print,'data lc_xspec.pha'
  print,'ign '
  print,where(tab.rate eq 0)+1
  print,'setpl ene'
  print,'cpd /xs'
  print,'pl ldat'
  print,'mo pow'
  print,'mo bknpow'
  print,'fit'
  print,'pl ldat res'
  print,'ipl'
  print,'line on 2'
  print,'pl'
  print,'quit'
end
