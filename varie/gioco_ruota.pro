pro gioco_ruota

                                ;Simulate game with 3 possible
                                ;outcomes, win 5 and continue, wins 10
                                ;and continues or game over


                                ;Simulate 1 game


  tot_win_array = lonarr(100000)

  for i = 1l, 100000 do begin
  
     tot_win = 0l
     game_over = 0l
     repeat begin
        a = fix(randomu(seed, 1)*3)
        if a eq 0 then tot_win = tot_win + 5
        if a eq 1 then tot_win = tot_win +10
        if a eq 2 then game_over = 1
     endrep until (game_over eq 1)
     tot_win_array[i-1] = tot_win
     
  endfor

  print, mean(tot_win_array)
  print, median(tot_win_array)
  plothist, tot_win_array
  
end
  
