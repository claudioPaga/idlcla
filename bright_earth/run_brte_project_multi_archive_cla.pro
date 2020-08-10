;10/8/12
;Andrew Spiers
pro run_brte_project_multi_archive_cla, obsid_file, g32
;Runs the brte_project32_archive version of the bright earth analysis on
;multiple datasets
;finding the w3 event files from the swift archive

if N_PARAMS() ne 2 then begin
    print, 'WARNING!!!!!!'
    print, 'Please enter the obsids table file and the grade selection (0= all grades, 1=g32 only)'
    print, "Ex: IDL> run_brte_project_multi32_archive, 'grbs_obsids.txt',1"
    print, '!!!!!!!!!!!!!'
    return
endif


;Reads in list of obs_ids
readcol, obsid_file, obsid, format = '(a)'

;Finds events, mkfs and sat files that match the obs_ids
;Runs brte_project32

for i= 0, n_elements(obsid)-1 do begin
    file_evt = FILE_SEARCH('/swift/data/archive/obs/'+obsid[i]+'/xrt/event/sw'+obsid[i]+'xpcw3po_uf.evt.gz',count = count_evt)
    mkf_match = file_search('/swift/data/archive/obs/'+obsid[i]+'/auxil/sw'+obsid[i]+'s.mkf.gz', count = count_mkf)
    coodin_match =  file_search('/swift/data/archive/obs/'+obsid[i]+'/auxil/sw'+obsid[i]+'sat.fits.gz', count =count_coord)
    print,file_evt[0]
    print,mkf_match[0]
    print, coodin_match[0]
    print,''
    if strtrim(string(count_evt),2) eq '1' and strtrim(string(count_mkf),2) eq '1' and strtrim(string(count_coord),2) eq '1' then brte_project_cla, file_evt[0], mkf_match[0], coodin_match[0], g32
                                ; kkk = get_kbrd(10)
endfor
end
