function inflation, initial_value, rate, years
  
  total = initial_value
  for i = 1, years do begin
     total = total+total*rate/100.
     print, i, total
  endfor

  return, total
end

  
