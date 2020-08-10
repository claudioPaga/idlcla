pro cleanflares,p
  read_lcfit,'lc_fit_out_idl_int9.dat',pname,p,perror
  file=findfile('lc_newout_phil.txt')
  if file[0] eq 'lc_newout_phil.txt' then begin
     lc=create_struct('time',0d,$
                   'tstart',0d,$
                   'tstop',0d,$
                   'src_rate',0d,$
                   'src_rate_err',0d,$
                   'tot_hard',0d,$
                   'tot_hard_err',0d,$
                   'exptime',0d,$
                   'src_counts',0d,$
                   'back_ctrate',0d,$
                   'det_sig',0d,$
                   'psf_corr',0d,$
                   'junk1',0d,$
                   'type',0d,$
                   'rate1',0d,$
                   'rate2',0d,$
                   'rate3',0d,$
                   'rate1_err',0d,$
                   'rate2_err',0d,$
                   'rate3_err',0d,$
                   'hard1',0d,$
                   'hard2',0d,$
                   'hard1_err',0d,$
                   'hard2_err',0d,$
                   'src_rate2',0d,$
                   'back_rate',0d,$
                   'back_area_corr',0d,$y
                   'pu_corr',0d,$
                   'psf_corr2',0d,$
                   'junk3',0d,$
                   'e1cts',0d,$
                   'e2cts',0d,$
                   'e3cts',0d,$
                   'junk4',0d,$
                   'junk5',0d,$
                   'junk6',0d,$
                   'tot_back_cts',0d)
   
     n=numlines(file)
     
     lc=replicate(lc,n-2)
     openr,lun,file,/get_lun
     line=readline(lun)
     line=readline(lun)
     for i=0,n-3 do begin
        line=readline(lun,delim=' ')
        w=where(line ne '')
        line=line[w]
        if n_elements(line) gt 24 and line[0] ne 'time' then begin 
           line=double(line)
           
           lc[i].time=line[0]
           lc[i].tstart=line[1]
           lc[i].tstop=line[2]
           lc[i].src_rate=line[3]
           lc[i].src_rate_err=line[4]
           lc[i].tot_hard=line[5]
           lc[i].tot_hard_err=line[6]
           lc[i].exptime=line[7]
           lc[i].src_counts=line[8]
           lc[i].back_ctrate=line[9]
           lc[i].det_sig=line[10]
           lc[i].psf_corr=line[11]
           lc[i].junk1=line[12]
           lc[i].type=line[13]
           lc[i].rate1=line[14]
           lc[i].rate2=line[15]
           lc[i].rate3=line[16]
           lc[i].rate1_err=line[17]
           lc[i].rate2_err=line[18]
           lc[i].rate3_err=line[19]
           lc[i].hard1=line[20]
           lc[i].hard2=line[21]
           lc[i].hard1_err=line[22]
           lc[i].hard2_err=line[23]
           lc[i].src_rate2=line[24]
           lc[i].back_rate=line[25]
           lc[i].back_area_corr=line[26]
           lc[i].pu_corr=line[27]
           lc[i].psf_corr2=line[28]
           lc[i].junk3=line[29]
           lc[i].e1cts=line[30]
           lc[i].e2cts=line[31]
           lc[i].e3cts=line[32]
           lc[i].junk4=line[33]
           lc[i].junk5=line[34]
           lc[i].tot_back_cts=line[35]
        endif 
     endfor 
     close,lun
     free_lun,lun
     
     w=where(lc.time ne 0)
     lc=lc[w]
     
     
     ;At this point lc now has all the info from the phil lc file
     
     ;Flares exclusion
     
     filegti=findfile('flares_gtis.dat')
     if filegti eq 'flares_gtis.dat' then begin
        readcol,filegti,f1,f2,format='(d,d)'
        for i=0,n_elements(f1)-1 do begin
           w=where(lc.time lt f1[i] or lc.time gt f2[i])
           lc=lc[w]
        endfor 
     endif
     
     ;Initial steep and jet-break phases exclusion (if they exist)
     
     if n_elements(p) gt 4 then begin
        
        if p[1] gt 1. then begin ;Here I assume there is an initial steep decay
                                ;I cut everything before the end of the steep decay
           index=where(lc.time gt p[2])
           lc=lc[index]
        endif
     
        if (p[n_elements(p)-1] gt p[n_elements(p)-3]+0.5 and p[n_elements(p)-3] gt 1.) then begin
                                ;I cut everything after the jet break
           index=where(lc.time lt p[n_elements(p)-2])
           lc=lc[index]
        endif
     endif
     
     ;I write out the new filetered LC
        
     file='lc_newout_noflares.txt'
     openw,lun,file,/get_lun
     s=sort(lc.time)
     sp='          '
     printf,lun,'time    t_start    t_stop   src_rate   src_rate_err   tot_hard   tot_hard_err   exptime    src_counts  back_ctrate   det_sig   psf_corr   junk   junk   rate1   rate2   rate3  rate1_err   rate2_err    rate3_err   hard1   hard2   hard1_err   hard2_err src_rate   back_rate   back_area_corr   pu_corr    psf_corr  somejunk  e1_cts  e2_cts  e3_cts  somejunk somejunk somejunk tot_back_cts  hard_ratio_1  hard_ratio_2  hard_ratio_err1  hard_ratio_err2 mean_time sigma_time'
     printf,lun,'junk2'
     for i=0,n_elements(lc)-1 do begin 
        printf,lun,ntostr(lc[s[i]].time)+sp+ntostr(lc[s[i]].tstart)+sp+ntostr(lc[s[i]].tstop)+sp+ntostr(lc[s[i]].src_rate)+sp+ntostr(lc[s[i]].src_rate_err)+sp+ntostr(1.)+sp+ntostr(.1)+sp+ntostr(lc[s[i]].exptime)+sp+ntostr(lc[s[i]].src_counts)+sp+ntostr(0)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(-999.)+sp+ntostr(lc[s[i]].type)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(lc[s[i]].pu_corr)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(lc[s[i]].tot_back_cts)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)+sp+ntostr(0.)
        
     endfor 
     close,lun
     free_lun,lun
  endif else print,'WARNING: There is no cleaned file yet!!!!!!!!!'
end
