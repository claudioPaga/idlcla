function figurine, n_fig_total, n_in_packet, check_unique = check_unique

                                ;Simulate how many figurine I need to
                                ;complete an album
  
  album = intarr(n_fig_total)
  
  repeat begin

     if check_unique then begin
;Check that there are no duplicates in
                                ;a packet (what happens in reality)

        repeat begin
           n_figurine = floor(randomu(seed, n_in_packet, /uniform)*n_fig_total)+1


        endrep until (n_elements(uniq(n_figurine)) eq n_in_packet)

     endif else n_figurine = floor(randomu(seed, n_in_packet, /uniform)*n_fig_total)+1
     
        
     for i = 0, n_in_packet-1 do album[n_figurine[i]-1] = album[n_figurine[i]-1]+1

     index0 = where(album eq 0, n0)

   endrep until n0 eq 0

  ;print, 'N tot figurine = ', total(album)
  
  array_return = [total(album), max(album)]
;  print, 'N tot pacchetti = ', total(album)/n_in_packet
  return, array_return
end


function figurine_scambi, n_friends, n_fig_total, n_in_packet, check_unique = check_unique

                                ;Simulate how many figurine I need to
                                ;complete an album
  
  album = intarr(n_fig_total)
  
  repeat begin

     if check_unique then begin
;Check that there are no duplicates in
                                ;a packet (what happens in reality)

        repeat begin
           n_figurine = floor(randomu(seed, n_in_packet, /uniform)*n_fig_total)+1


        endrep until (n_elements(uniq(n_figurine)) eq n_in_packet)

     endif else n_figurine = floor(randomu(seed, n_in_packet, /uniform)*n_fig_total)+1
     
                                ;I sum the stickers until I have at
                                ;least equal to n_friends stickers for
                                ;all stickers (equilvalent of having
                                ;completed n_friends albums   
     for i = 0, n_in_packet-1 do  album[n_figurine[i]-1] = album[n_figurine[i]-1]+1
     index_friends = where(album le n_friends, le_n_friends)

  endrep until le_n_friends eq 0

;  print, 'N tot figurine = ', total(album)
  
  array_return = [total(album), max(album)]
;  print, 'N tot pacchetti = ', total(album)/n_in_packet
  return, array_return
end




pro simulate_figurine, runs, n_friends, n_fig_total, n_in_packet, check_unique = check_unique

  DEVICE,DECOMPOSED=0.
tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1

  stats_figurine = intarr(runs, 2)
  
  for i = 0, runs-1 do begin

     stats_figurine[i, *] = figurine(n_fig_total, n_in_packet, check_unique = check_unique)

  endfor


  plothist, stats_figurine[*,0]

  print,'Stats on total number of stickers'
  print, 'Mean=', mean(stats_figurine[*,0])
  print, 'Median=', median(stats_figurine[*,0])
  print, 'Min=', min(stats_figurine[*,0])
  print, 'Max=', max(stats_figurine[*,0])
  print, 'Expected from maths = ', n_fig_total*alog(n_fig_total)


  print,'Stats on doubles'
  print, 'Mean=', mean(stats_figurine[*,1])
  print, 'Median=', median(stats_figurine[*,1])
  print, 'Min=', min(stats_figurine[*,1])
  print, 'Max=', max(stats_figurine[*,1])


  entry_device = !d.name
  plotfilenameps = 'figurine_histo.ps'
  set_plot,'ps'
  device, filename = plotfilenameps, /color
  plothist, stats_figurine[*,0], color=0, background=1
  plothist, stats_figurine[*,0], color=3, /overplot
  
  device,/close  
  set_plot,entry_device


stop

                                ;Run simulations for friends

  
for i = 0, runs-1 do begin

     stats_figurine[i, *] = figurine_scambi(n_friends, n_fig_total, n_in_packet, check_unique = check_unique)

  endfor


  plothist, stats_figurine[*,0]

  print,'Stats on total number of stickers'
  print, 'Mean=', mean(stats_figurine[*,0])
  print, 'Median=', median(stats_figurine[*,0])
  print, 'Min=', min(stats_figurine[*,0])
  print, 'Max=', max(stats_figurine[*,0])


  print,'Stats on doubles'
  print, 'Mean=', mean(stats_figurine[*,1])
  print, 'Median=', median(stats_figurine[*,1])
  print, 'Min=', min(stats_figurine[*,1])
  print, 'Max=', max(stats_figurine[*,1])


  entry_device = !d.name
  plotfilenameps = 'figurine_friends_histo.ps'
  set_plot,'ps'
  device, filename = plotfilenameps, /color
  plothist, stats_figurine[*,0], color=0, background=1
  plothist, stats_figurine[*,0], color=3, /overplot
  
  device,/close  
  set_plot,entry_device





  

end





     

  
