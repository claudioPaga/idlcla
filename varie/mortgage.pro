pro mortgage, prestito=prestito, rate=rate, years=years, inflazione = inflazione
  
  if not keyword_set(prestito) then prestito = 100000.
  if not keyword_set(rate) then rate = 1.
  if not keyword_set(years) then years = 25
  if not keyword_set(inflazione) then inflazione = 1.
  
  tot_n_payments = 12. * years
  rate_monthly = rate*1./1200.
  
  monthly_payment = rate_monthly * prestito * (1.+rate_monthly)^tot_n_payments / [(1.+rate_monthly)^tot_n_payments-1.]
  print, 'Monthly payments: ', monthly_payment
  yearly_payment = monthly_payment*12.
  print, 'Yearly payments: ', yearly_payment
  total_interest = monthly_payment * tot_n_payments - prestito
  print, 'Total interest = ', total_interest

                                ;Interests calculated as above are
                                ;relative to how much the money is
                                ;worth today, but that value should
                                ;really be adjusted for inflation I
                                ;think as for example something that
                                ;costs £500 today will cost £505 next
                                ;year with a 1% yearly inflation

  inflazione_monthly = inflazione*1./1200.
  monthly_inflation = inflazione_monthly * prestito * (1.+inflazione_monthly)^tot_n_payments / [(1.+inflazione_monthly)^tot_n_payments-1.]
  payment_less_inflation = (monthly_payment-monthly_inflation)*12.*years
  print, 'Total interest adjusted for inflation (metodo 1) ', payment_less_inflation
  
  
  adjusted_yearly_value = yearly_payment
  total = 0.
  for i = 1, years do begin
     adjusted_yearly_value = adjusted_yearly_value* (1. - inflazione/100.)
;     print, i, adjusted_yearly_value 
     total = total + adjusted_yearly_value
;     print, i, adjusted_yearly_value, total, format='(i, d12.2, d12.2)'
  endfor
  repayment_adjusted = total
  print, 'Total interest adjusted for inflation (metodo 2)', repayment_adjusted - prestito

;Calculate it as a loss in value
  yearly_payment_ad = yearly_payment
  total = 0.
  for i = 1, years do begin
     yearly_payment_ad = yearly_payment_ad + yearly_payment_ad*inflazione/100.
     total = total + (yearly_payment_ad - yearly_payment)
     ;print, i, yearly_payment_ad, total, Format='(i, F10.2, F10.2)'

  endfor
  print, 'Total interest adjusted for inflation (metodo 3)', total_interest-total



  print, 'Example of inflation, something that costs £1000 now will cost this much in future years, with a 2% infla:'
  now = 1000.
  infla = 0.028
  for i = 1, 25 do begin
     now = now*(1+infla)
     print, i, now
  endfor
  

  

                                ;Derive payments adjusted for
                                ;inflation (that is, how much the
                                ;yearly payment would be worth
                                ;in today's money)
;    yearly_payment_adjusted = fltarr(years-1)
  
  
;  yearly_payment_adjusted[0] = yearly_payment - yearly_payment*rate
;  for i = 2, years do begin
;     ;correct this!!!!
;     yearly_payment_adjusted[i-1] = yearly_payment_adjusted[i-1] - yearly_payment*rate
     
;  endfor



  
  ;Derive how much il prestito is worth in years time because of inflation
;  total = prestito
;  for i = 1, years do begin
;     total = total+total*inflazione/100.
;     print, i, total
;  endfor
;  prestito_inflated = total
;  print, 'Prestito ', prestito_inflated
  

end




  
