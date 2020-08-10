pro group_spectra_list

g0files = file_search('*g0.pi')
g02files = file_search('*g0-2.pi')

for k = 0,n_elements(g0files)-1 do begin
    comando = '/swift/home/apb/work/python_scripts/fireball/grppha.py -m 1 -b _bgd '+g0files[k]+' '+g02files[k]
    spawn,comando
endfor
end
