pro badpix
readcol,'/bulk/yankees2/xrt/pass_data/pass_20051170649/pass_20051170649_science.0_LDP20283.pc.cat',pc_name,format='(a)'

hot_pix_detx=[0]
hot_pix_dety=[0]

totale_eventi=0L
for ii=0,n_elements(pc_name)-1 do begin
    tab=mrdfits('/bulk/yankees2/xrt/pass_data/pass_20051170649/'+pc_name[ii],1,hdr)
    hot_pix_detx=[hot_pix_detx,tab.detx]
    hot_pix_dety=[hot_pix_dety,tab.dety]
    n_eventi=total(tab.detx)
    totale_eventi=totale_eventi+n_eventi
endfor
;print,hot_pix_detx
print,'Totale eventi: ',totale_eventi

print,'List of hot pixels:'

contar_rep=0
for i=0,150 do begin
    conta_rep=0
    for iii=0,n_elements(hot_pix_detx)-1 do begin
        if ((hot_pix_detx[iii] eq hot_pix_detx[i]) and (hot_pix_dety[iii] eq hot_pix_dety[i])) then conta_rep=conta_rep+1
    endfor
    if (conta_rep ge 80) then begin
        print,hot_pix_detx[i],hot_pix_dety[i]
     ;   print,''
    endif
endfor

print,'Numero di PC frmes: ',n_elements(pc_name)
print,'Numero di eventi: ',n_elements(hot_pix_detx),n_elements(hot_pix_dety)

end

