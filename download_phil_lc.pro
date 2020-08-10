pro download_phil_lc
  readcol,'listadir.txt',dir,obs,format='(a,a)'
  dir='GRB'+dir
  obs='00'+obs  
;  if n_params() eq 0 then begin
;     print,'syntax - download_phil_lc,dir'
;     print,"         (e.g. dir='GRB080123')"
;     return
;  endif 
  
;  cd,!mdata
;  if n_elements(dir) eq 0 then dir=file_search('GRB*')
  for i=0,n_elements(dir)-1 do begin 
     scritta='mkdir '+dir(i)
     spawn,scritta
     cd,dir[i]
    ; spawn,'rm WTCURVE.qdp'
    ; spawn,'rm PCCURVE.qdp'
     url='http://www.swift.ac.uk/xrt_curves/'+obs[i]+'/WTCURVE.qdp'
     com='wget '+url
     print,com
     spawn,com
     url='http://www.swift.ac.uk/xrt_curves/'+obs[i]+'/PCCURVE.qdp'
     com='wget '+url
     print,com
     spawn,com
     
     cd,'..'
  endfor
  
  return
end 
