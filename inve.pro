pro inve

read,inve_start,prompt='Investimento iniziale: '
read,inve_month,prompt='Investimento mensile (primi 3 anni): '
read,tasso,prompt='Tasso guadagno fondo (percentuale medio): '
read,inve_money,prompt='Investimento Money market mensile (primi 3 anni): '
totale_cA=inve_start*(1.-0.0575)
mensile_cA=inve_month*(1-0.0575)
totale_cB=inve_start
mensile_cB=inve_month
money_tot=0.

print,''

for i=1,10 do begin
    if (i lt 4) then begin  ;Add money montly for first 3 years 
        tassoanno=0.
        for m=1,12 do begin
            tassomensile=(((randomn(seed,1)*5)*0.01)+(tasso*0.01))/12.0 ;faccio fluttuare il tasso mensile
                                ;attorno al tasso che dico io di piu' o meno il 10%
                                ;il tasso che dico io pero' lo divido
                                ;per 12, perche' se ho un tasso del
                                ;12% annuo vuol dire che e' dell'1% mensile
                                ;Introdotto accumulo mensile 
            totale_cA=(totale_cA+mensile_cA)*(1+tassomensile)
            totale_cB=(totale_cB+mensile_cB)*(1+tassomensile)
            tassoanno=tassoanno+tassomensile
        endfor
        totanno_cA=totale_cA*(1-0.0077) ;0.0077=tasse per classe A
        totanno_cB=totale_cB*(1-0.0154)
        totaleinvestito=inve_start+(12*inve_month)*i+inve_money*12*i
        if (i eq 1) then totalepagato=totanno_cB-inve_start*((6-i)*0.01)-(12*inve_month)*5*0.01 
        if (i eq 2) then totalepagato=totanno_cB-inve_start*((6-i)*0.01)-(12*inve_month)*(5+4)*0.01
        if (i eq 3) then totalepagato=totanno_cB-inve_start*((6-i)*0.01)-(12*inve_month)*(5+4+3)*0.01
    endif else begin
        tassoanno=(randomn(seed,1)*5)*0.01+tasso*0.01
        ;tassoanno=tasso*0.01
        totale_cA=totale_cA*(1+tassoanno)
        totale_cB=totale_cB*(1+tassoanno)
        totanno_cA=totale_cA*(1-0.0077) ;0.0077=tasse per classe A
        totanno_cB=totale_cB*(1-0.0154) ;0.0154=tasse per classe B
        if (i eq 4) then totalepagato=totanno_cB-inve_start*((6-i)*0.01)-(12*inve_month)*(4+3+2)*0.01
        if (i eq 5) then totalepagato=totanno_cB-inve_start*((6-i)*0.01)-(12*inve_month)*(3+2+1)*0.01
        if (i eq 6) then totalepagato=totanno_cB-(12*inve_month)*(2+1)*0.01
        if (i eq 7) then totalepagato=totanno_cB-(12*inve_month)*0.01
        if (i ge 8) then totalepagato=totanno_cB
    endelse
    if i lt 4 then begin
        money_tot=(money_tot+(inve_money*12))*1.04
    endif else money_tot=money_tot*1.04

    pr1=string(strtrim(tassoanno,2))+' , '+string(strtrim(totaleinvestito,2))
    print,'Inflazione/Totale investito anno ',strtrim(i,2),': ',pr1

    conto_cA=totanno_cA+money_tot
    conto_cB=totalepagato+money_tot
    pr2=string(strtrim(conto_cA,2))+' , '+string(strtrim(conto_cB,2))
    print,'Totale fondoA/fondoB: ',pr2
    print,''
endfor


inflazione=totaleinvestito
infla=0.03
print,'Inflazione del 3%: '
for k=1,10 do begin
    inflazione=inflazione*(1.03)
    print,'Anno ',k,' Infla: ',inflazione
endfor
print,'Prezzo dopo 10 anni dell investimento: ',strtrim(inflazione,2)

end

         
