pro be_copy_files,inobs

readcol,inobs,obsids,format='(a)'

for i = 0,n_elements(obsids)-1 do begin
    strobs = strtrim(obsids[i],2)
    stringcom = 'cp /swift/data/archive/obs/'+strobs+'/xrt/event/sw'+strobs+'xwtw2po_cl.evt.gz .'
    spawn,stringcom
    stringcom = 'cp /swift/data/archive/obs/'+strobs+'/auxil/sw'+strobs+'sao.fits.gz .'
    spawn,stringcom
    stringcom = 'cp /swift/data/archive/obs/'+strobs+'/auxil/sw'+strobs+'s.mkf.gz .'
    spawn,stringcom
endfor
spawn,'gunzip *gz'

listeve = file_search('sw*evt')
listsao = file_search('sw*sao.fits')
listmkf = file_search('sw*mkf')

for j = 0, n_elements(listeve)-1 do begin
    plot_radiator_earth_angle_satfile,listeve[j],listsao[j],listmkf[j]
    print,'Next? (type to continue)'
    k = get_kbrd(10)
endfor

end
