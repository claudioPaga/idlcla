function cti_correction_function, e_measured, ac, ctife, ctibeta, cmax, threshold

  cti = ctife * (e_measured)^(ctibeta)

  counter = 0

  e_iter = e_measured

  losses = 0.
  
  repeat begin

     losses_previous = losses
     losses = e_iter * (1. - (1.-cti)^ac)
     e_iter = e_measured + losses
     counter = counter + 1
     cti = ctife * (e_iter)^(ctibeta)
     
     ;print, 'Counter, losses, e_iter, cti'
     ;print, counter, losses, e_iter, cti
     
  endrep until (counter eq cmax) or (losses - losses_previous lt threshold)


  return, e_iter

end
