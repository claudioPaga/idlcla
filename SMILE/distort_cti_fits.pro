FUNCTION event_splitting, events_list_energies_keV

  ;;; Input - the list of original X-ray energies, in KeV
  ;;; Returns - the list of split X-ray energies, in KeV, returned as
  ;;;           an array of [9,n_elements(events_list)]

  ;;; Splits the energies of each event based on Georges' table
  ;;; of single/doubles/triples/quadruples probabilities

  ;;; Grade splitting info

   split_energy_ev = [100,150.0,200.0,250.0,300.0,350.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,1250.0,1500.0,1750.0,1800.0,1838,1840.0,1900.0,1950.0,2000.0,2500.0,3000.0,3500.0,4000.0,4500.0,5000.0,5500.0,6000.0]
   singles = [1.0,0.99944496,0.99246104	,0.972028602	,0.951572566	,0.935876794	,0.922449224	,0.90012	,0.88703	,0.88074	,0.880257605	,0.882645967	,0.889813599	,0.904895411	,0.910234142	,0.911368865	,0.913263473	,0.820366407	,0.817104997	,0.816808102	,0.817316251	,0.818198182	,0.83563004	,0.856692606	,0.869273067	,0.87486931	,0.881303336	,0.882872356	,0.881611937	,0.881277892]
   doubles =  [0.0	,0.00055504	,0.00753896	,0.027971398	,0.04841743	,0.064083204	,0.077360774	,0.09916	,0.11158	,0.11725	,0.117392348	,0.114673255	,0.1075092	,0.092571941	,0.087334121	,0.086247985	,0.084134719	,0.172153443	,0.175035404	,0.175366264	,0.17484273	,0.173791738	,0.156767762	,0.136603495	,0.124458375	,0.119324961	,0.112935666	,0.110901534	,0.112925253	,0.113173243]
   triples = [0.0, 0.0, 0.0, 0.0, 1.00E-05, 4.00E-05, 0.000160002, 0.0006	,0.0011, 0.00151, 0.001650033, 0.001820528, 0.001644423, 0.001459139, 0.001505361, 0.001515298 ,0.001581209 ,0.004190084 ,0.004531661, 0.004613322	,0.004630602	, 0.004730047	,0.004142095	,0.003849866	,0.003211261	,0.003296474	,0.00314628	,0.003071325	,0.00287409	,0.002786239]
   quadruples = [0.0, 0.0, 0.0, 0.0,	0.0, 0.0, 3.00E-05 ,0.00012	,0.00029	,0.0005	,0.000700014	,0.000860249	,0.001032778	,0.001073509	,0.000926376	,0.000867852	,0.001020599	,0.003290066	,0.003327938	,0.003212313	,0.003210417	,0.003280033	,0.003460104	,0.002854034	,0.003057297	,0.002509256	,0.002614717	,0.003154785	,0.00258872	,0.002762627]

   ;;; Array to store the energies in the 9 pixles X-ray event 3x3
  ;;; pixel window
  ;;; Top left pixel index = 0, bottom right pixel index = 8
  ;;; Pixels 0, 1, 2 won't contribute to the downstream charge as it would
  ;;; be outside the window
  totEvents = n_elements(events_list_energies_keV)
  list_energy_grid = dblarr(9, totEvents)
  
  for countEvents = 0, totEvents - 1 do begin
     ;;; Determine if the event is a single, double, triple or
     ;;; quadruple drawing from a random uniform distrib comparing
     ;;; with George's grade probability tables.
     
     ran = randomu(seed, 1)
     resultSingle = INTERPOL(singles, split_energy_ev, events_list_energies_keV[countEvents]*1000.)
     resultDoubles = INTERPOL(doubles, split_energy_ev, events_list_energies_keV[countEvents]*1000.)
     resultTriples = INTERPOL(triples, split_energy_ev, events_list_energies_keV[countEvents]*1000.)
     ;print, ran, events_list_energies_keV[countEvents]*1000., resultSingle, resultSingle+ resultDoubles, resultSingle+ resultDoubles+ resultTriples
     ;stop
     CASE 1 OF
        ;;; SINGLE
        (ran le  resultSingle): list_energy_grid[4,countEvents] = events_list_energies_keV[countEvents]

        ;;; DOUBLE
        (ran le  resultSingle+resultDoubles): BEGIN
           splitEn1 = events_list_energies_keV[countEvents] * randomu(seed, 1)
           splitEn2 = events_list_energies_keV[countEvents] - splitEn1
           list_energy_grid[4, countEvents] = max([splitEn1, splitEn2])
           position = fix(4.*randomu(seed, 1))
           CASE position OF
              0: list_energy_grid[1,countEvents] = min([splitEn1, splitEn2]) ; Double up
              1: list_energy_grid[5,countEvents] = min([splitEn1, splitEn2]) ; Double right
              2: list_energy_grid[7,countEvents] = min([splitEn1, splitEn2]) ; Double down
              3: list_energy_grid[3,countEvents] = min([splitEn1, splitEn2]) ; Double left
           ENDCASE
        END
        
        ;;;TRIPLE
        (ran le  resultSingle+resultDoubles+resultTriples): BEGIN
           ;;; Pick the maximum by ensuring its energy is more than
           ;;; half the total energy
           central = 0.5*events_list_energies_keV[countEvents] + 0.5*events_list_energies_keV[countEvents]*randomu(seed, 1)
           list_energy_grid[4, countEvents] = central
           
           splitEn1 = (events_list_energies_keV[countEvents] - central) * randomu(seed, 1)
           splitEn2 = events_list_energies_keV[countEvents] - central - splitEn1

           ;;; Pick the pattern at random by the possible 4 valid
           ;;; 3-pixel orientations 
           position = fix(4.*randomu(seed, 1))
           CASE position OF
              0: list_energy_grid[[1,5], [countEvents,countEvents]] = [splitEn1, spliten2] ; Up and right
              1: list_energy_grid[[5,7], [countEvents,countEvents]] = [splitEn1, spliten2] ; Right and Down
              2: list_energy_grid[[7,3], [countEvents,countEvents]] = [splitEn1, spliten2] ; Down and left
              3: list_energy_grid[[3,1], [countEvents,countEvents]] = [splitEn1, spliten2] ; Left and Up
           END
        END

        ;;; QUADRUPLE
        ELSE: BEGIN
           ;;; Pick the maximum by ensuring its energy is more than
           ;;; half the total energy
           central = 0.5*events_list_energies_keV[countEvents] + 0.5*events_list_energies_keV[countEvents]*randomu(seed, 1)
           list_energy_grid[4, countEvents] = central
           ;;; Distribute the remaining energy within the other 3
           ;;; pixels
           second = (events_list_energies_keV[countEvents] - central) * randomu(seed, 1)
           third =  (events_list_energies_keV[countEvents] - central - second) *randomu(seed, 1)
           fourth = events_list_energies_keV[countEvents] - central  - second - third

           position = fix(4.*randomu(seed, 1))

           CASE position OF
              0: list_energy_grid[[1,2,5], [countEvents,countEvents,countEvents]] = [second, third, fourth] ; Top right box
              1: list_energy_grid[[5,7,8], [countEvents,countEvents,countEvents]] = [second, third, fourth] ; Bottom right box
              2: list_energy_grid[[3,6,7], [countEvents,countEvents,countEvents]] = [second, third, fourth] ; Bottom left box
              3: list_energy_grid[[0,1,3], [countEvents,countEvents,countEvents]] = [second, third, fourth] ; Top left obx
           END
        END
     ENDCASE
  endfor

  return, list_energy_grid

end


FUNCTION event_splitting_smart, events_list_energies_keV, split_ll_keV

  ;;; Input - events_list_energies_keV = the list of original X-ray
  ;;;         energies, in KeV
  ;;;         split_ll = the energy lower limit for a valid energy in
  ;;;         a pixel (lower is considered noise)
  ;;; Returns - the list of split X-ray energies, in KeV, returned as
  ;;;           an array of [9,n_elements(events_list)]
  ;;; Note - try a more efficient way to split the events compared to FUNCTION event_splitting
  
  ;;; Splits the energies of each event based on Georges' table
  ;;; of single/doubles/triples/quadruples probabilities

  ;;; Grade splitting info
   split_energy_ev = [100,150.0,200.0,250.0,300.0,350.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,1250.0,1500.0,1750.0,1800.0,1838,1840.0,1900.0,1950.0,2000.0,2500.0,3000.0,3500.0,4000.0,4500.0,5000.0,5500.0,6000.0]
   singles = [1.0,0.99944496,0.99246104	,0.972028602	,0.951572566	,0.935876794	,0.922449224	,0.90012	,0.88703	,0.88074	,0.880257605	,0.882645967	,0.889813599	,0.904895411	,0.910234142	,0.911368865	,0.913263473	,0.820366407	,0.817104997	,0.816808102	,0.817316251	,0.818198182	,0.83563004	,0.856692606	,0.869273067	,0.87486931	,0.881303336	,0.882872356	,0.881611937	,0.881277892]
   doubles =  [0.0	,0.00055504	,0.00753896	,0.027971398	,0.04841743	,0.064083204	,0.077360774	,0.09916	,0.11158	,0.11725	,0.117392348	,0.114673255	,0.1075092	,0.092571941	,0.087334121	,0.086247985	,0.084134719	,0.172153443	,0.175035404	,0.175366264	,0.17484273	,0.173791738	,0.156767762	,0.136603495	,0.124458375	,0.119324961	,0.112935666	,0.110901534	,0.112925253	,0.113173243]
   triples = [0.0, 0.0, 0.0, 0.0, 1.00E-05, 4.00E-05, 0.000160002, 0.0006	,0.0011, 0.00151, 0.001650033, 0.001820528, 0.001644423, 0.001459139, 0.001505361, 0.001515298 ,0.001581209 ,0.004190084 ,0.004531661, 0.004613322	,0.004630602	, 0.004730047	,0.004142095	,0.003849866	,0.003211261	,0.003296474	,0.00314628	,0.003071325	,0.00287409	,0.002786239]
   quadruples = [0.0, 0.0, 0.0, 0.0,	0.0, 0.0, 3.00E-05 ,0.00012	,0.00029	,0.0005	,0.000700014	,0.000860249	,0.001032778	,0.001073509	,0.000926376	,0.000867852	,0.001020599	,0.003290066	,0.003327938	,0.003212313	,0.003210417	,0.003280033	,0.003460104	,0.002854034	,0.003057297	,0.002509256	,0.002614717	,0.003154785	,0.00258872	,0.002762627]

   ;;; Array to count number of splits
   countSplitsArray = [0,0,0]
   
   ;;; Array to store the energies in the 9 pixles X-ray event 3x3
  ;;; pixel window
  ;;; Top left pixel index = 0, bottom right pixel index = 8
  ;;; Pixels 0, 1, 2 won't contribute to the downstream charge as it would
  ;;; be outside the window
  totEvents = n_elements(events_list_energies_keV)
  list_energy_grid = dblarr(9, totEvents)
  
  for countEvents = 0, totEvents - 1 do begin
     ;;; Determine if the event is a single, double, triple or
     ;;; quadruple drawing from a random uniform distrib comparing
     ;;; with George's grade probability tables.
     
     ran = randomu(seed, 1)    
     resultSingle = INTERPOL(singles, split_energy_ev, events_list_energies_keV[countEvents]*1000.)
     resultDoubles = INTERPOL(doubles, split_energy_ev, events_list_energies_keV[countEvents]*1000.)
     resultTriples = INTERPOL(triples, split_energy_ev, events_list_energies_keV[countEvents]*1000.)
     
     CASE 1 OF
        ;;; SINGLE
        (ran le  resultSingle): list_energy_grid[4,countEvents] = events_list_energies_keV[countEvents]
   
        ;;; DOUBLE
        (ran le  resultSingle+resultDoubles): BEGIN
           ;;; Check if enough energy in PHA for a split, if so, split
           ;;; energy between pixels assuring its above the split
           ;;; threshold split_ll_keV
           over_split = events_list_energies_keV[countEvents] - 2.*split_ll_keV
           if over_split gt 0. then begin
              countSplitsArray[0] +=1
              splitEn = randomu(seed,2)
              splitEn = splitEn/total(splitEn) * over_split + split_ll_keV
              list_energy_grid[4, countEvents] = max(splitEn, indexmax)
              position = fix(4.*randomu(seed, 1))
              CASE position OF
                 0: list_energy_grid[1,countEvents] = min(splitEn) ; Double up
                 1: list_energy_grid[5,countEvents] = min(splitEn) ; Double right
                 2: list_energy_grid[7,countEvents] = min(splitEn) ; Double down
                 3: list_energy_grid[3,countEvents] = min(splitEn) ; Double left
              ENDCASE
           endif else list_energy_grid[4, countEvents] = events_list_energies_keV[countEvents] ; If not enough energy for a split keep it a single event.
        END
        
        ;;;TRIPLE
        (ran le  resultSingle+resultDoubles+resultTriples): BEGIN
           ;;; Check if enough energy in PHA for a split, if so, split
           ;;; energy between pixels assuring its above the split
           ;;; threshold split_ll_keV
           over_split = events_list_energies_keV[countEvents] - 3.*split_ll_keV
           if over_split gt 0. then begin
              countSplitsArray[1] +=1
              splitEn = randomu(seed,3)
              splitEn = splitEn/total(splitEn) * over_split + split_ll_keV         
              centralValue = max(splitEn, indexMax)
              list_energy_grid[4, countEvents] = centralValue
           
           ;;; Indices of split pixels (bit convoluted, possibly there
           ;;; is a better way to get these indices.)         
              indexNoMax = [0 ne indexMax, 1 ne indexMax, 2 ne indexMax]
              indexNoMax = where(indexNoMax eq 1)
           
           ;;; Pick the pattern at random by the possible 4 valid
           ;;; 3-pixel orientations 
           position = fix(4.*randomu(seed, 1))
           CASE position OF
              0: list_energy_grid[[1,5], [countEvents,countEvents]] = splitEn[indexNoMax] ; Up and right
              1: list_energy_grid[[5,7], [countEvents,countEvents]] = splitEn[indexNoMax] ; Right and Down
              2: list_energy_grid[[7,3], [countEvents,countEvents]] = splitEn[indexNoMax] ; Down and left
              3: list_energy_grid[[3,1], [countEvents,countEvents]] = splitEn[indexNoMax] ; Left and Up
           END
           endif else list_energy_grid[4, countEvents] = events_list_energies_keV[countEvents] ; If not enough energy for a split keep it a single event.
        END

        ;;; QUADRUPLE
        ELSE: BEGIN
           ;;; Check if enough energy in PHA for a split, if so, split
           ;;; energy between pixels assuring its above the split
           ;;; threshold split_ll_keV
           over_split = events_list_energies_keV[countEvents] - 4.*split_ll_keV
           if over_split gt 0. then begin
              countSplitsArray[2] +=1
              splitEn = randomu(seed,4)
              splitEn = splitEn/total(splitEn) * over_split + split_ll_keV          
              centralValue = max(splitEn, indexMax)
              list_energy_grid[4, countEvents] = centralValue
           
           ;;; Indices of split pixels (bit convoluted, possibly there
           ;;; is a better way to get these indices.)         
              indexNoMax = [0 ne indexMax, 1 ne indexMax, 2 ne indexMax, 3 ne indexMax]
              indexNoMax = where(indexNoMax eq 1)

              ;;; Pick the pattern at random (4 valid 4-pixels orientations)
              position = fix(4.*randomu(seed, 1))
              CASE position OF
                 0: list_energy_grid[[1,2,5], [countEvents,countEvents,countEvents]] = splitEn[indexNoMax] ; Top right box
                 1: list_energy_grid[[5,7,8], [countEvents,countEvents,countEvents]] = splitEn[indexNoMax] ; Bottom right box
                 2: list_energy_grid[[3,6,7], [countEvents,countEvents,countEvents]] = splitEn[indexNoMax] ; Bottom left box
                 3: list_energy_grid[[0,1,3], [countEvents,countEvents,countEvents]] = splitEn[indexNoMax] ; Top left obx
              END
              endif else list_energy_grid[4, countEvents] = events_list_energies_keV[countEvents]; If not enough energy for a split keep it a single event.
           END
     ENDCASE
  endfor

  print, 'Singles, Doubles, Triples, Quadruples = ' , totEvents - total(countSplitsArray), ' ', countSplitsArray[0], ' ', countSplitsArray[1],' ', countSplitsArray[2]
  return, list_energy_grid

end






function cti_process, events_list_energies_keV, rawx, rawy, f_parallel_trap_density_pixel, f_serial_trap_density_pixel, f_cti_alpha, f_cti_beta, f_serial_cti_alpha, f_serial_cti_beta, f_key_release

;;;
;;; Input = List of energies in a specific pixel position of the 3x3
;;; window
;;; Output = 3 x n_events array that includes:
;;;          Damaged pixels
;;;          Parallel downstream pixel trailing charge
;;;          Serial downstream pixel trailing charge
;;;
 ; COMMON shareCTI, trap_species_fraction_parallel, trap_species_release_constant_parallel, transfer_time_parallel, trap_species_fraction_serial, trap_species_release_constant_serial, transfer_time_serial, cti_alpha, cti_beta, cti_corr_alpha, cti_corr_beta, serial_cti_alpha, serial_cti_beta, serial_cti_corr_alpha, serial_cti_corr_beta, pixel_binning, key_release

  ;;; Initialise trap parameters
  pixel_binning = 6
  trap_species_fraction_parallel = [0.6, 0.4]
  trap_species_release_constant_parallel = [1., 10.]
  transfer_time_parallel = 1.
  trap_species_fraction_serial = [0.5, 0.1]
  trap_species_release_constant_serial = [1., 10.]
  transfer_time_serial = 1.

  
  list_energies = events_list_energies_keV
    
  xraysTot = n_elements(list_energies)
  list_energies_cti = dblarr(xraysTot)
  losses_array_ev = dblarr(xraysTot)
  parallel_recovered_kev_stored = dblarr(xraysTot)
  parallel_downstream_kev = dblarr(xraysTot)
  parallel_downstream_kev_stored = dblarr(xraysTot)
  serial_recovered_kev_stored = dblarr(xraysTot)
  serial_downstream_kev = dblarr(xraysTot)
  serial_downstream_kev_stored = dblarr(xraysTot)
  serial_cti_damaged_stored = dblarr(xraysTot)
  
 ;;; *** PARALLEL TRANSFER, IMAGE FRAME'

  
 ;;; Simple model of electron capture model.
 ;;; capture_prob = prob of capturing an electron from a charge packet
 ;;; of size list_energies
 ;;; losses_mean = expected electron losses, given the capture
 ;;; probability, the number of transfers and the density of traps in
 ;;; a pixel
 ;;; losses_sigma = standard deviation of expected caputes (model the
 ;;; probabilistic nature of trapping)
 
 capture_prob = 1. - exp(-f_cti_alpha * list_energies^(1-f_cti_beta))

                                ;Distribution of number of transfers
                                ;(this corresponds to the DETY
                                ;position of the X-ray event. The
                                ;assumption is of a uniform
                                ;distribution along the column.
 ;;; Transform the input DETY from binned Row Coordinates to Pixel
 ;;; coordinates, randomising its location within the 6 pixels, and
 ;;; assuming all charge is in a pixel (not realistic)
 
 n_transfers = rawy*6 + fix(randomu(Seed, xraysTot)*6) + 1
 ;;; Losses in electrons
 losses_mean = n_transfers * capture_prob * f_parallel_trap_density_pixel 
 losses_sigma = sqrt(n_transfers * f_parallel_trap_density_pixel * capture_prob * (1. - capture_prob))
 random_array = randomn(seed, n_elements(list_energies))
 ;;; Losses in eV
 losses_array_ev = (random_array * losses_sigma + losses_mean) * 3.65 

                                ;Avoid "negative" losses (can be
                                ;caused by low losses_mean value and a
                                ;large losses_sigma value
 index_l0 = where(losses_array_ev lt 0., nl0)
 if nl0 gt 0 then losses_array_ev[index_l0] = 0
 
;;; Apply parallel losses, avoid negative energies replacing negative
;;; values with zeros and replacing the losses with the total energy
;;; in that pixel
 list_energies_cti  = (list_energies*1000. - losses_array_ev)/1000. ; Damaged flux in keV
 index_neg = where(list_energies_cti lt 0., nneg)
 if nneg gt 0 then begin
    list_energies_cti[index_neg] = 0.0
    losses_array_ev[index_neg] = list_energies[index_neg]*1000.
 endif
 list_energies_cti_pre_recovery = list_energies_cti 
 
  if f_key_release then begin
    ;;; Number of pixels [0-5] within which released charge will be
    ;;; within the 6 pixel binning and added back to the initial energy.
    pixels_released = n_transfers mod pixel_binning
    parallel_en_losses_kev = list_energies - list_energies_cti
                                ;Sum recoveries over the trap species
    for counter_species_parallel = 0, n_elements(trap_species_fraction_parallel) - 1 do begin
       parallel_charge_released_within_binning =  trap_species_fraction_parallel[counter_species_parallel]  * parallel_en_losses_kev * (1. - exp(-pixels_released*transfer_time_parallel/ trap_species_release_constant_parallel[counter_species_parallel]))
       parallel_recovered_kev_stored += parallel_charge_released_within_binning
    endfor
    list_energies_cti  += parallel_recovered_kev_stored
 endif
 
    ;;; ********** Release of charge in downstream 6x-binned pixel **********

    ;;; To contribute to downstream pixel, the released charge will
    ;;; need to be released within [pixel_released - pixel_released+6]
    ;;; Formula (1-exp(-pixel_released + 6 /tau)) - (1-exp(-pixel_released/tau)) 
    ;;; = exp(-pixel_released/ tau) - exp(-pixel_released+6/ tau)

    ;;; Sum over the trap species of the charge released in downstream pixels
    for counter_species_parallel = 0, n_elements(trap_species_fraction_parallel) - 1 do begin
       parallel_charge_released_downstream =  trap_species_fraction_parallel[counter_species_parallel]  * parallel_en_losses_kev * (exp(-pixels_released* transfer_time_parallel/ trap_species_release_constant_parallel[counter_species_parallel]) - exp(-(6+pixels_released)* transfer_time_parallel/ trap_species_release_constant_parallel[counter_species_parallel]))
       parallel_downstream_kev += parallel_charge_released_downstream
    endfor
    
    ;;; Charge at this stage is transferred to the store section ,
    ;;; accumulated over 6 pixels and transferred over 719 more rows.

    ;;; list_energies_cti already represents the energy in the 6
    ;;; pixels and simply needs to be trasferred to the bottom of the
    ;;; store section

    ;;; *********** Apply store section parallel CTI losses ********

    capture_prob = 1. - exp(-f_cti_alpha * list_energies_cti^(1-f_cti_beta))

                                ;Distribution of number of transfers
                                ;(this corresponds to the DETY
                                ;position of the X-ray event. The
                                ;assumption is of a uniform
                                ;distribution along the column.
 ;;; Transform the input DETY from binned Row Coordinates to Pixel
 ;;; coordinates, randomising its location within the 6 pixels, and
 ;;; assuming all charge is in a pixel (not realistic)
 
 n_transfers_framestore = 719
 ;;; Losses in electrons
 losses_mean = n_transfers_framestore * capture_prob * f_parallel_trap_density_pixel 
 losses_sigma = sqrt(n_transfers_framestore * f_parallel_trap_density_pixel * capture_prob * (1. - capture_prob))
 random_array = randomn(seed, n_elements(list_energies_cti))
 ;;; Losses in eV
 losses_array_ev = (random_array * losses_sigma + losses_mean) * 3.65 

                                ;Avoid "negative" losses (can be
                                ;caused by low losses_mean value and a
                                ;large losses_sigma value
 index_l0 = where(losses_array_ev lt 0., nl0)
 if nl0 gt 0 then losses_array_ev[index_l0] = 0
 
 ;;; Apply parallel losses, avoid negative energies
 
 list_energies_cti_framestore  = (list_energies_cti*1000. - losses_array_ev)/1000. ; Damaged flux in keV
 index_neg = where(list_energies_cti_framestore lt 0., nneg)
 if nneg gt 0 then begin
    list_energies_cti_framestore[index_neg] = 0.0
    losses_array_ev[index_neg] = list_energies_cti[index_neg]*1000.
 endif
 list_energies_cti_pre_recovery_framestore = list_energies_cti_framestore

  ;;; Losses (keV) during frame store parallel transfers.
  parallel_en_losses_framestore_kev = list_energies_cti - list_energies_cti_framestore

   
  ;;; Downstream pixel - Derive electrons released in downstream pixel
  ;;;                    from charge lost during store parallel
  ;;;                    transfers and add it to the flux in the
  ;;;                    downstream pixel.   

  ;;; Sum over the trap species of the charge released in downstream
  ;;; pixel
  ;;; In store section this should simply correspond to the following 1
  ;;; pixel, so it's charge released within 1 transfer
  for counter_species_parallel = 0, n_elements(trap_species_fraction_parallel) - 1 do begin
     parallel_charge_released_downstream =  trap_species_fraction_parallel[counter_species_parallel]  * parallel_en_losses_framestore_kev * (1. - exp(-transfer_time_parallel/ trap_species_release_constant_parallel[counter_species_parallel]))
     parallel_downstream_kev += parallel_charge_released_downstream
  endfor

  ;;; list_energies_cti_framestore ==> energies after image + store sections
  ;;; parallel_downstream_kev  ==> charge accumulated in downstream pixel.
 
 ;;; Apply serial CTI distortion

 ;;; Charge is transferred using two nodes, so a max of 2250 serial
 ;;; transfers occour

 ;;; Transform the input DETY from binned Row Coordinates to Pixel
 ;;; coordinates, randomising its location within the 6 pixels, and
 ;;; assuming all charge is in a pixel (not realistic)
 ;;; Calculate the number of serial transfers from the closest node,
 ;;; store it in serial_n_transfers.

  RAWX_unbinned = rawx * 6 + fix(randomu(Seed, xraysTot)*6) + 1
  serial_n_transfers = intarr(xraysTot)
  for xrayCounter = 0, xraysTot-1 do serial_n_transfers[xrayCounter] = min([RAWX_unbinned[xrayCounter], 4510-RAWX_unbinned[xrayCounter]])


 serial_capture_prob = 1. - exp(-f_serial_cti_alpha * list_energies_cti_framestore^(1-f_serial_cti_beta))
 serial_losses_mean = serial_n_transfers * serial_capture_prob * f_serial_trap_density_pixel 
 serial_losses_sigma = sqrt(serial_n_transfers * f_serial_trap_density_pixel * serial_capture_prob * (1. - serial_capture_prob))
 serial_random_array = randomn(seed, n_elements(list_energies_cti))
 serial_losses_array_ev = (serial_random_array * serial_losses_sigma + serial_losses_mean) * 3.65
                                  ;Avoid "negative" losses (can be
                                ;caused by low losses_mean value and a
                                ;large losses_sigma value
 index_l0 = where(serial_losses_array_ev lt 0., nl0)
 if nl0 gt 0 then serial_losses_array_ev[index_l0] = 0
 
 ;;; Apply serial losses, avoid negative energies
 list_energies_cti  = (list_energies_cti_framestore*1000. - serial_losses_array_ev)/1000.
 index_neg = where(list_energies_cti lt 0., nneg)
 if nneg gt 0 then begin
    list_energies_cti[index_neg] = 0.0
    serial_losses_array_ev[index_neg] = list_energies_cti[index_neg]*1000.
 endif
 
 list_energies_cti_pre_serial_recovery = list_energies_cti_framestore
 
; Recovery of release charge in binning
 if f_key_release then begin
    pixels_released = serial_n_transfers mod pixel_binning
    serial_en_losses_kev = list_energies_cti_framestore - list_energies_cti
                                ;Sum recoveries over the trap species
    for counter_species_serial = 0, n_elements(trap_species_fraction_serial) - 1 do begin
       serial_charge_released_within_binning =  trap_species_fraction_serial[counter_species_serial]  * serial_en_losses_kev * (1. - exp(-pixels_released*transfer_time_serial/ trap_species_release_constant_serial[counter_species_serial]))
       serial_recovered_kev_stored += serial_charge_released_within_binning
    endfor
    list_energies_cti  += serial_recovered_kev_stored
 endif
 
  ;;; ********** Release of charge in downstream serial 6x-binned pixel **********

    ;;; To contribute to downstream pixel, the released charge will
    ;;; need to be released within [pixel_released - pixel_released+6]
    ;;; Formula (1-exp(-pixel_released + 6 /tau)) - (1-exp(-pixel_released/tau)) 
    ;;; = exp(-pixel_released/ tau) - exp(-pixel_released+6/ tau)

    ;;; Sum over the trap species of the charge released in downstream pixels
 for counter_species_serial = 0, n_elements(trap_species_fraction_serial) - 1 do begin
    serial_charge_released_binning = trap_species_fraction_serial[counter_species_serial]  * serial_en_losses_kev * (exp(-pixels_released*transfer_time_serial/ trap_species_release_constant_serial[counter_species_serial]) - exp(-(6+pixels_released)*transfer_time_serial/ trap_species_release_constant_serial[counter_species_serial]))
    serial_downstream_kev += serial_charge_released_binning
 endfor
    
;;; serial_downstream_kev is the charge in the binned downstream
;;; serial pixel (binned pixel to the right/left of central pixel,
;;; depending on the node).
 

 ;;; Derive serial losses of downstream electrons (electrons in the
 ;;; binned pixel just above)

 serial_capture_prob = 1. - exp(-f_serial_cti_alpha * parallel_downstream_kev^(1-f_serial_cti_beta))
 serial_losses_mean = serial_n_transfers * serial_capture_prob * f_serial_trap_density_pixel 
 serial_losses_sigma = sqrt(serial_n_transfers * f_serial_trap_density_pixel * serial_capture_prob * (1. - serial_capture_prob))
 serial_random_array = randomn(seed, n_elements(parallel_downstream_kev))
 serial_losses_array_ev = (serial_random_array * serial_losses_sigma + serial_losses_mean) * 3.65
                                  ;Avoid "negative" losses (can be
                                ;caused by low losses_mean value and a
                                ;large losses_sigma value
 index_l0 = where(serial_losses_array_ev lt 0., nl0)
 if nl0 gt 0 then serial_losses_array_ev[index_l0] = 0
 
 ;;; Apply serial losses, avoid negative energies 
 serial_downstream_damaged_kev = (parallel_downstream_kev*1000.- serial_losses_array_ev)/1000.
 findNegIndex = where(serial_downstream_damaged_kev lt 0, negTot)
 if negTot gt 0 then begin
    serial_downstream_damaged_kev[findNegIndex] = 0.0
    serial_losses_array_ev[index_neg] = parallel_downstream_kev[index_neg]*1000.
 endif
  
 ;;; Recovery of release charge in binning downstream pixel
 if f_key_release then begin
    pixels_released = serial_n_transfers mod pixel_binning
    serial_en_losses_kev = parallel_downstream_kev - serial_downstream_damaged_kev
                                ;Sum recoveries over the trap species
    serial_recovered_kev_stored = dblarr(xraysTot)
    for counter_species_serial = 0, n_elements(trap_species_fraction_serial) - 1 do begin
       serial_charge_released_within_binning =  trap_species_fraction_serial[counter_species_serial]  * serial_en_losses_kev * (1. - exp(-pixels_released*transfer_time_serial/ trap_species_release_constant_serial[counter_species_serial]))
       serial_recovered_kev_stored += serial_charge_released_within_binning
    endfor
    serial_downstream_damaged_kev  += serial_recovered_kev_stored
 endif
 ;findNegIndex = where(serial_downstream_damaged_kev lt 0, negTot)
 ;if negTot gt 0. then  serial_downstream_damaged_kev[findNegIndex] = 0.0;
  
 ;;; serial_downstream_damaged_kev = the charge in the downstream
 ;;; pixel (above the central pixel) after serial losses + release in
 ;;; the binned pixel

 
 ;;; Remove X-rays with energy in central pixel below event threshold
 ;;; (zero to start with) due to losses 

  ;;; Replace negative energies if present with 0 energy
  index_negative_energies = where(list_energies_cti lt 0., n_negative)
  if n_negative gt 0 then begin
     list_energies_cti[index_negative_energies] = 0.
     serial_downstream_kev[index_negative_energies] = 0.
     serial_downstream_damaged_kev[index_negative_energies] = 0. 
  endif
  serial_cti_damaged_stored = list_energies_cti

  ; Report if CTI too high and all energies lost.
 index_positive_energies = where(list_energies_cti gt 0., n_positive)
 ;print, 'N events, N positive after serial readout: ', n_elements(list_energies_cti), n_positive
 ;if n_positive eq 0 then begin
 ;   print, 'CTI very high, all events are lost!!!'
    ;return
 ;endif
 
 ;;; *** REPORT ON EVENT ENERGIES AFTER PARALLEL+SERIAL CTI DISTORTION ****

 ;;; Damaged energies at this stage:
 ;;; serial_cti_damaged_stored ==> Damaged flux in central pixel,
 ;;; binned, after parallel image+store section damage, serial damage,
 ;;; and recovery from released charge within binned pixel from both
 ;;; parallel and serial transfers
 ;;; serial_downstream_kev ==> Charge released in binned pixel after
 ;;; central pixel.
 ;;; serial_downstream_damaged_kev ==> Charge released in downstream binned pixel
 ;;; during parallel transfer, serial-cti- damaged with electrons
 ;;; recovered from charge
 ;;; released within binned
 ;;; pixel during serial transfer.

 ;;; The combination of these three charges is the binned split event.
 ;;; serial_cti_damaged_stored == Central pixel
 ;;; serial_downstream_kev == Trailing serial pixel
 ;;; serial_downstream_damaged_kev == Trailing parallel pixel

 damaged_array = dblarr(3,xraysTot)
 damaged_array[0,*] = serial_cti_damaged_stored
 damaged_array[1,*] = serial_downstream_damaged_kev
 damaged_array[2,*] = serial_downstream_kev
                                ; damaged_array = [serial_cti_damaged_stored, serial_downstream_damaged_kev, serial_downstream_kev]
  
  return, damaged_array

end



function cti_correction, input_damaged_list

  ;;; !!!!!!!  NOTE - JUST A PLACEHOLDER FOR NOW !!!!!!!!!!!
  ;;;
  ;;; Code from previous version of the code, need fixing!!!!


  return, input_damaged_list


 ;;; Attempt the derivation of improved CTI correction coefficients
 ;;; based on simulated damages
 
 !P.Multi = [0, 1, 2, 0, 0]

  cti_parallel = 1. - (parallel_cti_damaged_stored/list_energies_stored)^(1./n_transfers_stored)
  cti_serial = 1. - (serial_cti_damaged_stored/serial_list_energies)^(1./serial_n_transfers_stored)

                                ; Parallel CTI linear fit
                                ; Avoid CTI = 0 in fit, re-evaluate
                                ; CTI-correction parameters only if
                                ; there are enough events for a good
                                ; fit, otherwise use the input ones.
  index_values_ok_parallel = where(parallel_n_transfers_stored gt 100. and cti_parallel ne 0, n_fit_values_parallel)

  if n_fit_values_parallel gt 1000 then begin
  
     cti_parallel_logscale = alog10(cti_parallel[index_values_ok_parallel])
     list_energies_logscale_ev = alog10(list_energies_stored[index_values_ok_parallel]*1000.)
     res_linear_parallel = linfit(list_energies_logscale_ev, cti_parallel_logscale, yfit = yfit_linfit_parallel)
     print, 'Parallel CTI fit: ', res_linear_parallel
     plot, list_energies_logscale_ev, cti_parallel_logscale, psym=1, color=0, background=1, ytit = 'Parallel CTI' 
     oplot, list_energies_logscale_ev, yfit_linfit_parallel, color=3

     print, 'Warning, new parameters for parallel CTI corrections were derived'
     print, 'Old values = ', cti_corr_alpha, cti_corr_beta
     cti_corr_alpha = 10^res_linear_parallel[0]
     cti_corr_beta = res_linear_parallel[1]
     print, 'New values = ', cti_corr_alpha, cti_corr_beta
         
  endif
     
     
                                ; Serial CTI linear fit
                                ; Avoid CTI = 0 in fit, re-evaluate
                                ; CTI-correction parameters only if
                                ; there are enough events for a good
                                ; fit, otherwise use the input ones.
  index_values_ok_serial = where(serial_n_transfers_stored gt 100. and cti_serial ne 0, n_fit_values_serial)
  if n_fit_values_serial gt 1000 then begin
   
     parallel_cti_damaged_stored_logscale_ev = alog10(serial_list_energies[index_values_ok_serial]*1000.)
     cti_serial_logscale = alog10(cti_serial[index_values_ok_serial])
     res_linear_serial = linfit(parallel_cti_damaged_stored_logscale_ev, cti_serial_logscale, yfit = yfit_linfit_serial)
     print, 'Serial CTI fit: ', res_linear_serial
     plot, parallel_cti_damaged_stored_logscale_ev, cti_serial_logscale, psym=1, color=0, background=1, ytit = 'Serial CTI'
     oplot, parallel_cti_damaged_stored_logscale_ev, yfit_linfit_serial, color=3

     print, 'Warning, new parameters for Serial CTI corrections were derived'
     print, 'Old values = ', serial_cti_corr_alpha, serial_cti_corr_beta
     serial_cti_corr_alpha = 10^res_linear_serial[0]
     serial_cti_corr_beta = res_linear_serial[1]
     print, 'New values = ', serial_cti_corr_alpha, serial_cti_corr_beta
                 
  endif
  
!P.Multi = [0, 1, 1, 0, 0]

    
; Remove events below the lower level discriminator, these events will be lost
 index_detection_ll_energies = where(list_energies_cti gt detection_ll/1000., n_ll) ;detection_ll is in eV, but list_energies_cti is in keV!
 print, 'N events, N above threshold after serial losses, total lost events: ', n_elements(list_energies_cti), n_ll,  n_elements(list_energies)-n_ll
 if n_ll gt 0 then begin
    serial_cti_damaged_stored_ll = list_energies_cti[index_detection_ll_energies]
    serial_n_transfers_stored_ll = serial_n_transfers[index_detection_ll_energies]
    parallel_n_transfers_stored_ll = parallel_n_transfers_stored[index_detection_ll_energies]
    list_energies_ll = list_energies_stored[index_detection_ll_energies]
 endif else begin
    print, 'CTI very high, all events are lost below the low level discriminator!!!'
 endelse


 ; CTI Corrections
 ; A loop is needed as the corection function works on single energy events.
                                ; Estimate event location (due to 6pixels binning)
 list_energies_cti_corrected = dblarr(n_elements(serial_cti_damaged_stored_ll))
 cti_corrections_ev_storage = dblarr(n_elements(serial_cti_damaged_stored_ll))
 for count_events = 0, n_elements(list_energies_cti_corrected) - 1 do begin
    ; Serial correction first
    eve_energy_ev = serial_cti_damaged_stored_ll[count_events] * 1000.
     ; Estimate event location (due to 6pixels binning)
    x_location = floor((floor(serial_n_transfers_stored_ll[count_events]/pixel_binning)+randomu(seed, 1))*pixel_binning+1.)
    e_corr_serial_cti = cti_correction_function(eve_energy_ev, x_location, serial_cti_corr_alpha, serial_cti_corr_beta, cmax, threshold)

 ; Parallel correction (takes binning
 ; into account as when correcting we
 ; only have 6-pixel binning info   
    y_location = floor((floor(parallel_n_transfers_stored_ll[count_events]/pixel_binning)+randomu(seed, 1))*pixel_binning+1.)
    e_corr_parallel_cti = cti_correction_function(e_corr_serial_cti, y_location, cti_corr_alpha, cti_corr_beta, cmax, threshold)
    list_energies_cti_corrected[count_events] = e_corr_parallel_cti/1000.
    cti_corrections_ev_storage[count_events] = e_corr_parallel_cti - eve_energy_ev
    ;if losses_array_ev[count_events]-cti_corrections_ev_storage[count_events] gt 30. then stop
 endfor

 plothist, list_energies, xhist_input, yhist_input, bin=hw[0]*2, background=1, color=0, boxplot=0
 plothist, list_energies_cti, xhist, yhist, bin=hw[0]*2, /overplot, color=3, boxplot=0
 plothist, list_energies_cti_corrected, xhist_corrected, yhist_corrected, bin=hw[0]*2, /overplot, color=4, boxplot=0


                                ;Determine ymax for postscipt plot
 ymaxplot = max([yhist_input, yhist, yhist_corrected])
 plothist, list_energies, xhist_input, yhist_input, bin=hw[0]*2, background=1, color=0, yr=[0, ymaxplot], boxplot=0
 plothist, serial_cti_damaged_stored_ll, xhist, yhist, bin=hw[0]*2, /overplot, color=3, boxplot=0
 plothist, list_energies_cti_corrected, xhist_corrected, yhist_corrected, bin=hw[0]*2, /overplot, color=4, boxplot=0
 string_split = strsplit(qdp_filename, '.qdp', /regex, /extract)

                                ;Input CTI parameters are arbitrary,
                                ;if set too high all charge can be
                                ;lost. CTI qdp file and plot only
                                ;generated if there is some positive
                                ;chage (at least 10 bins)
 only_positive = where(xhist ge 0, n_positive)

 if n_positive gt 10 then begin

    xhist_gt0 = xhist[only_positive]
    yhist_gt0 = yhist[only_positive]

                                ;Estimate the errors in the
                                ;cti-distorted fluxes, based on the
                                ;best-fit parameters of the input
                                ;errors vs counts values
    err_cti_kev = sqrt(res_lin[0] + res_lin[1] * yhist_gt0)
 
    outfilenameqdp = 'cti_'+string_split[0]+'_detxmin'+strtrim(string(detx_min[0]),2)+'_detxmax'+strtrim(string(detx_max[0]),2)+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'_release'+strtrim(string(key_release),2)+'.qdp'
    openw, lu, outfilenameqdp, /get_lun
    print_array = [transpose(xhist_gt0), transpose(replicate(hw[0],n_elements(xhist_gt0))), transpose(yhist_gt0/(hw[0] * 2 * expo_qdp[0])), transpose(err_cti_kev)]
    printf, lu, print_array
    free_lun, lu

    ;Save cti-corrected spectrum
    xhist_corrected_gt0 = xhist_corrected[only_positive]
    yhist_corrected_gt0 = yhist_corrected[only_positive]
    err_cti_corrected_kev = sqrt(res_lin[0] + res_lin[1] * yhist_corrected_gt0)
    
    outfilenameqdp = 'cti_with_correction_'+string_split[0]+'_detxmin'+strtrim(string(detx_min[0]),2)+'_detxmax'+strtrim(string(detx_max[0]),2)+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'_release'+strtrim(string(key_release),2)+'.qdp'
    openw, lu, outfilenameqdp, /get_lun
    print_array = [transpose(xhist_corrected_gt0), transpose(replicate(hw[0],n_elements(xhist_corrected_gt0))), transpose(yhist_corrected_gt0/(hw[0] * 2 * expo_qdp[0])), transpose(err_cti_corrected_kev)]
    printf, lu, print_array
    free_lun, lu
    


    
    if keyword_set(key_plot) then begin
       entry_device = !d.name
       plotfilenameps = 'cti_with_correction_'+string_split[0]+'_detxmin'+strtrim(string(detx_min[0]),2)+'_detxmax'+strtrim(string(detx_max[0]),2)+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'_release'+strtrim(string(key_release),2)+'.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color
       titstring = 'Input [black] + distorted [red] + corrected [blue] '+ string_split[0]
       plothist, list_energies, xhist, yhist, bin=hw[0]*2, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "Cnts/s/keV", yr=[0, ymaxplot], boxplot=0, charsize=0.8,xr=[0.05,3]
       plothist, serial_cti_damaged_stored_ll, bin=hw[0]*2, /overplot, color=3, boxplot=0
       plothist, list_energies_cti_corrected, bin=hw[0]*2, /overplot, color=4, boxplot=0
       expostr = 'Expo = '+ strtrim(string(floor(expo_qdp[0])), 2)+ ' s'
       xyouts, 0.7, 0.8, expostr, /norm, charsize=0.8
       driftstr = 'Drift = '+ strtrim(string(floor(dety_max[0]-dety_min[0])), 2)+ ' pixels'
       xyouts, 0.7, 0.75, driftstr, /norm, charsize=0.8
       xtrmaxstr = 'DETX = '+ strtrim(string(floor(detx_min[0])), 2)+ '-' +  strtrim(string(floor(detx_max[0])), 2)
       xyouts, 0.7, 0.7, xtrmaxstr, /norm, charsize=0.8
       ytrmaxstr = 'DETY = '+ strtrim(string(floor(dety_min[0])), 2)+ '-' +  strtrim(string(floor(dety_max[0])), 2)
       xyouts, 0.7, 0.65, ytrmaxstr, /norm, charsize=0.8
       device,/close  
       set_plot,entry_device

       
       plotfilenameps = 'cti_with_correction_'+string_split[0]+'_detxmin'+strtrim(string(detx_min[0]),2)+'_detxmax'+strtrim(string(detx_max[0]),2)+'_detymin'+strtrim(string(dety_min[0]),2)+'_detymax'+strtrim(string(dety_max[0]),2)+'_release'+strtrim(string(key_release),2)+'_log.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color
       titstring = 'Input [black] + distorted [red] + corrected [blue] ' + string_split[0]
       plothist, list_energies, bin=hw[0]*2, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "Cnts/s/keV",xr=[0.05,3],/ylog, yr=[1, ymaxplot], xst=1, yst=1, boxplot=0, charsize=0.8
       plothist, serial_cti_damaged_stored_ll, bin=hw[0]*2, /overplot, color=3, boxplot=0
       plothist, list_energies_cti_corrected, bin=hw[0]*2, /overplot, color=4, boxplot=0
       xyouts, 0.7, 0.8, expostr, /norm, charsize=0.8
       xyouts, 0.7, 0.75, driftstr, /norm, charsize=0.8
       xyouts, 0.7, 0.7, xtrmaxstr, /norm , charsize=0.8
       xyouts, 0.7, 0.65, ytrmaxstr, /norm, charsize=0.8
       device,/close  
       set_plot,entry_device

       ;Additional plots of distribution of losses vs DETY and ENERGY

       !P.Multi = [0, 1, 2, 0, 0]
        
       plotfilenameps = 'losses_vs_dety_energy_'+string_split[0]+'.ps'
       set_plot,'ps'
       device, filename = plotfilenameps, /color
       titstring = 'CTI losses '+ string_split[0]
       plot, list_energies, losses_array_ev, background=1, color=0, tit = titstring, xtit="Energy [keV]", ytit = "CTI Losses (eV)",xr=[0.05,3], psym=1
       plot, n_transfers, losses_array_ev, background=1, color=0, tit = titstring, xtit="DETY", ytit = "CTI Losses (eV)", psym=1
       
       device,/close  
       set_plot,entry_device


        !P.Multi = [0, 1, 1, 0, 0]
       
       
    endif

 endif else print, 'WARNING - CTI TOO HIGH, ALL CHARGE IS LOST!!!!!'

  

  return, list_energies_cti_corrected

end



pro distort_cti_fits, fits_filename, expo, parallel_trap_density_pixel, serial_trap_density_pixel, cti_alpha, cti_beta, cti_corr_alpha, cti_corr_beta, serial_cti_alpha, serial_cti_beta, serial_cti_corr_alpha, serial_cti_corr_beta, cmax, threshold, key_release = key_release, key_plot = key_plot, detection_ll_eV = detection_ll_eV, split = split, key_cti_corr = key_cti_corr

;;;
;;; CP, 19 Dec 2018
;;; Summary - Distorts Xray events from a fits file using a simple CTI model
;;;           Applies correction to distorted spectrum, using original
;;;           CTI parameters used to apply the distortion, that is,
;;;           CTI(Fe) = 10^-4 (SMILE EoL)
;;;           CTI beta = -0.25 (same as Swift)
;;;
;;;  
;;;  
;;; Input: fits_filename - X-ray events fits file
;;;        expo - Exposure time of fits file
;;;        parallel_trap_density_pixel = number of traps in a pixel in
;;;        image section
;;;        serial_trap_density_pixel = number of traps in a pixel in
;;;        the readout register 
;;;        cti_alpha, cti_beta = parallel CTI distortion parameters
;;;        cti_corr_alpha, cti_corr_beta = CTI correction pars, not the
;;;        same as distortion, from SMILE EoL
;;;        preditions,
;;;        serial_cti_alpha, serial_cti_beta = serial CTI distortion parameters
;;;        serial_cti_corr_alpha, serial_cti_corr_beta = CTI correction pars, not the
;;;        same as distortion, from SMILE EoL
;;;        preditions,  
;;;        cmax, threshold = CTI correction iteration parameters, cmax
;;;        = max number of trials in iteration, threshold = min
;;;        differenct between derived value of corrected energy in the
;;;        iteration, in eV  
;;;        split = If set, simulate split (non-isolated) events based
;;;        on George's info.
;;;  
;;; KEYWORD PARAMETERS:
;;;        key_release = recover energy within the 6-pixel binning,
;;;        default = 0
;;;        key_plot = generate plots
;;;        detection_ll_eV = Min eV value to detect events (like lower  
;;;        level discriminator for XRT, it's the value to
;;;        "separate" real X-ray events from noise. Default = 0 (no threshold)
;;;  
;;; Output: CTI distorted fits file
;;;  
;;; HISTORY - Written by CP, 14 Aug 2020
;;;  
;;;           Modified from distort_cti_qdp_dety_min_max
;;;           added procedure to correct for CTI losses after CTI
;;;           damage, to check how well energies can be recovered
;;;    
;;;           Modified from distort_cti_qdp_pixel.pro
;;;           Distortion over a range of dety values, dety_min, dety_max
;;;           Distortion assumes observations are uniformely
;;;           distributed between dety_min and dety_max
;;;
;;;          CP, 18/1/2019
;;;          Modified from distort_cti_qdp_dety_min_max to implement serialCTI
;;;          distortion for a source imaged between detx_min and detx_max
;;;          and by added low threshold discriminator (keyword detection_ll) to
;;;          exclude low-energy X-rays that the CCD will not be able to
;;;          detect. it's the low level discriminator in XRT, to
;;;          separate noise from real X-rays  
;;;
;;;          CP, 14 Aug 2020
;;;          Modified from distort_cti_qdp_detx_dety_min_max_and_cti_correction.pro 
;;;          Uses input fits file instead of qdp file
;;;          Option to generate non-isolated X-rays (non grade0)
;;;          Binning 6x6 implemented
;;;          Charge release treatment implemented
  
  
   if n_params() lt 14 then begin
    print,"Please enter correct parameters: distort_cti_fits, fits_filename, expo, trap_density_pixel_parallel, trap_density_pixel_serial, cti_alpha, cti_beta, cti_corr_alpha, cti_corr_beta, serial_cti_alpha, serial_cti_beta, serial_cti_corr_alpha, serial_cti_corr_beta, cmax, threshold"
    print,"Please enter correct parameters: distort_cti_fits, 'P100001SXI0033IMEVLI0000.FIT', 25200, 1., 1., 0.031, 0.44, 0.018843, -0.7, 0.012, 0.49, 0.00604860, -0.7, 10, 1, key_release = 1, key_plot=1, detection_ll_eV=100.0, split = 1"
    return
 endif  
  
  DEVICE,DECOMPOSED=0.
  tvlct,[255,0,255,0,0],[255,255,0,0,0],[255,0,0,255,0],1

  ;;; Initialise keywords if needed.
  
  if not(keyword_set(key_release)) then key_release = 0 ; Default = no recovery of released charge in pixel binning
  if not(keyword_set(split)) then split = 0 ; Default = no splits, only grade 0 events

  ;;; Output fits name
  splitString = strsplit(fits_filename, '.', /extract)
  fits_filename_output = splitString[0]+'CTI.FIT'
  ps_filename_output = splitString[0]+'CTI.ps'

  ;;; Read in fits file
  tab = mrdfits(fits_filename,1,hd1)
  tab0 = mrdfits(fits_filename,0,hd0)
  xraysTot = 0l
  xraysTot = n_elements(tab.PI)
  
  ;;; ***************   EVENT SPLITTING  ************
  ;;;
  ;;; Call the event_splitting function to split the energy of the
  ;;; events between the allowed X-ray grades.
  ;;; The returned array is a 9 x N_tot_events array, with the 9
  ;;; representing the pixels in the 3x3 window.
  ;;; If no grade splitting requested assign all energy to the central
  ;;; pixels


  if not(keyword_set(detection_ll_eV)) then split_ll_keV = 0 else split_ll_keV = detection_ll_eV/1000.

  ;;; Initialise array to include the 3x3 energies.
  list_energy_grid = dblarr(9, xraysTot)
  if split then list_energy_grid = event_splitting_smart(tab.PHA/1000., split_ll_keV) else list_energy_grid[4,*] = tab.PHA/1000.
  ;if split then list_energy_grid = event_splitting(tab.PHA/1000.) else list_energy_grid[4,*] = tab.PHA/1000.
  
  
  ;;; Array to accumulate charge of CTI processed pixels
  ;;; The 3by3 window has indices:
  ;;;
  ;;;   0 1 2
  ;;;   3 4 5
  ;;;   6 7 8
  ;;;
  ;;; So for example index 4 is relative to the central pixel, index 1
  ;;; to its parallel trailing charge 
  
  damaged_9pixels_array_ev = dblarr(9, xraysTot)

  ;;; Process pixels of position 0, accumulate charge in damaged array
  processed_pixels_array = cti_process(list_energy_grid[0,*], tab.RAWX, tab.RAWY, parallel_trap_density_pixel, serial_trap_density_pixel, cti_alpha, cti_beta, serial_cti_alpha, serial_cti_beta, key_release)
  damaged_9pixels_array_ev[0,*] = processed_pixels_array[0,*]*1000.
  damaged_9pixels_array_ev[1,*] = processed_pixels_array[2,*]*1000. ;Serial trailing charge added

  ;;; Process pixels of position 1, accumulate charge in damaged array
  processed_pixels_array = cti_process(list_energy_grid[1,*], tab.RAWX, tab.RAWY, parallel_trap_density_pixel, serial_trap_density_pixel, cti_alpha, cti_beta, serial_cti_alpha, serial_cti_beta, key_release)
  damaged_9pixels_array_ev[1,*] += processed_pixels_array[0,*]*1000.
  damaged_9pixels_array_ev[2,*] = processed_pixels_array[2,*]*1000. ;Serial trailing charge added
  
  ;;; Process pixels of position 2, accumulate charge in damaged array
  processed_pixels_array = cti_process(list_energy_grid[2,*], tab.RAWX, tab.RAWY, parallel_trap_density_pixel, serial_trap_density_pixel, cti_alpha, cti_beta, serial_cti_alpha, serial_cti_beta, key_release)
  damaged_9pixels_array_ev[2,*] += processed_pixels_array[0,*]*1000.

  ;;; Process pixels of position 3, accumulate charge in damaged array
  processed_pixels_array = cti_process(list_energy_grid[3,*], tab.RAWX, tab.RAWY, parallel_trap_density_pixel, serial_trap_density_pixel, cti_alpha, cti_beta, serial_cti_alpha, serial_cti_beta, key_release)
  damaged_9pixels_array_ev[3,*] = processed_pixels_array[0,*]*1000.
  damaged_9pixels_array_ev[0,*] += processed_pixels_array[1,*]*1000. ;Parallel trailing charge
  damaged_9pixels_array_ev[4,*] = processed_pixels_array[2,*]*1000. ;Serial trailing charge
  
  ;;; Process pixels of position 4, accumulate charge in damaged array
  processed_pixels_array = cti_process(list_energy_grid[4,*], tab.RAWX, tab.RAWY, parallel_trap_density_pixel, serial_trap_density_pixel, cti_alpha, cti_beta, serial_cti_alpha, serial_cti_beta, key_release)
  damaged_9pixels_array_ev[4,*] += processed_pixels_array[0,*]*1000.
  damaged_9pixels_array_ev[1,*] += processed_pixels_array[1,*]*1000. ;Parallel trailing charge
  damaged_9pixels_array_ev[5,*] = processed_pixels_array[2,*]*1000.  ;Serial trailing charge
  
  ;;; Process pixels of position 5, accumulate charge in damaged array
  processed_pixels_array = cti_process(list_energy_grid[5,*], tab.RAWX, tab.RAWY, parallel_trap_density_pixel, serial_trap_density_pixel, cti_alpha, cti_beta, serial_cti_alpha, serial_cti_beta, key_release)
  damaged_9pixels_array_ev[5,*] += processed_pixels_array[0,*]*1000.
  damaged_9pixels_array_ev[2,*] += processed_pixels_array[1,*]*1000. ;Parallel trailing charge

;;; Process pixels of position 6, accumulate charge in damaged array
  processed_pixels_array = cti_process(list_energy_grid[6,*], tab.RAWX, tab.RAWY, parallel_trap_density_pixel, serial_trap_density_pixel, cti_alpha, cti_beta, serial_cti_alpha, serial_cti_beta, key_release)
  damaged_9pixels_array_ev[6,*] = processed_pixels_array[0,*]*1000.
  damaged_9pixels_array_ev[3,*] += processed_pixels_array[1,*]*1000. ;Parallel trailing charge
  damaged_9pixels_array_ev[7,*] = processed_pixels_array[2,*]*1000.  ;Serial trailing charge

;;; Process pixels of position 7, accumulate charge in damaged array
  processed_pixels_array = cti_process(list_energy_grid[7,*], tab.RAWX, tab.RAWY, parallel_trap_density_pixel, serial_trap_density_pixel, cti_alpha, cti_beta, serial_cti_alpha, serial_cti_beta, key_release)
  damaged_9pixels_array_ev[7,*] += processed_pixels_array[0,*]*1000.
  damaged_9pixels_array_ev[4,*] += processed_pixels_array[1,*]*1000. ;Parallel trailing charge
  damaged_9pixels_array_ev[8,*] = processed_pixels_array[2,*]*1000.  ;Serial trailing charge
  
;;; Process pixels of position 8, accumulate charge in damaged array
  processed_pixels_array = cti_process(list_energy_grid[8,*], tab.RAWX, tab.RAWY, parallel_trap_density_pixel, serial_trap_density_pixel, cti_alpha, cti_beta, serial_cti_alpha, serial_cti_beta, key_release)
  damaged_9pixels_array_ev[8,*] += processed_pixels_array[0,*]*1000.
  damaged_9pixels_array_ev[5,*] += processed_pixels_array[1,*]*1000. ;Parallel trailing charge
  
  ;;; Low detection limit filtering
  ;;; NOTE - The event reconstruction shouldn't be part of the
  ;;;        CTI processing really. All charge should be reported in
  ;;;        the PHAS and filtered accordingly to the low level
  ;;;        distcriminator in the grade/event assignment.
  ;;;        But in case the low limit is set, for ease of checking
  ;;;        the effect of a low level limit the energies of the
  ;;;        pixels in the window are set to zero if below the limit,
  ;;;        apart from the central pixel that is always kept as
  ;;;        calculated after the CTI processing..

  if keyword_set(detection_ll_eV) then begin
     for pixel_index = 0, 8 do begin
        if pixel_index ne 4 then begin
           for events_counter = 0l, xraysTot - 1 do begin
              if damaged_9pixels_array_ev[pixel_index,events_counter] lt detection_ll_eV then  damaged_9pixels_array_ev[pixel_index,events_counter] = 0.
           endfor
        endif
     endfor
  endif
  
  damaged_total_array_ev = total(damaged_9pixels_array_ev, 1)

  
  mwrfits, tab0, fits_filename_output, hd0, /create
  
  newtab = {EVENTS, TIME: tab[0].TIME, RAWX: tab[0].RAWX, RAWY: tab[0].RAWY, DETX: tab[0].DETX, DETY: tab[0].DETY, X: tab[0].X, Y: tab[0].Y , PHA:tab[0].PHA, PI:tab[0].PI, PHAS:list_energy_grid[*,0]*1000., PHASCTI:damaged_9pixels_array_ev[*,0]}
  outtab = replicate({EVENTS}, xraysTot)
  for count= 0l, xraysTot-1 do outtab[count] = {TIME: tab[count].TIME, RAWX: tab[count].RAWX, RAWY: tab[count].RAWY, DETX: tab[count].DETX, DETY: tab[count].DETY, X: tab[count].X, Y: tab[count].Y ,PHA:tab[count].PHA, PI:tab[count].PI, PHAS:list_energy_grid[*,count]*1000., PHASCTI:damaged_9pixels_array_ev[*,count]}
  mwrfits, outtab, fits_filename_output, hd1

  plothist, tab.PHA, /xlog, /ylog
  plothist, damaged_total_array_ev, /overplot, color=2
  plothist, damaged_9pixels_array_ev[4,*], /overplot, color=3
  xyouts, 0.6, 0.85, 'Input PHA', /norm, color = 0, charsize = 1.2
  xyouts, 0.6, 0.8, 'Summed CTI-damaged PHAS', /norm, color = 2, charsize = 1.2
  xyouts, 0.6, 0.75, 'Central CTI-damaged PHA', /norm, color = 3, charsize = 1.2

  if keyword_set(key_plot) then begin
       entry_device = !d.name
       set_plot,'ps'
       device, filename = ps_filename_output, /color
       ;titstring = 'Input [black] + distorted [red] + corrected [blue] '+ string_split[0]
       plothist, tab.PHA, /xlog, /ylog
       plothist, damaged_total_array_ev, /overplot, color=2
       plothist, damaged_9pixels_array_ev[4,*], /overplot, color=3
       xyouts, 0.6, 0.85, 'Input PHA', /norm, color = 0
       xyouts, 0.6, 0.8, 'Summed CTI-damaged PHAS', /norm, color = 2
       xyouts, 0.6, 0.75, 'Central CTI-damaged PHA', /norm, color = 3
       device,/close  
       set_plot,entry_device

    endif


  ;;; CORRECT CTI DAMAGE

  if keyword_set(key_cti_corr) then cti_corrected_list = cti_correction(damaged_9pixels_array_ev)

  
  stop
 
 

 
 ;en6kev = 6.
 ;losses6keV= (1. - exp(-cti_alpha * en6kev^(1-cti_beta))) * n_transfers * trap_density_pixel * 3.65
 ;en1kev = 1.
 ;losses1keV= (1. - exp(-cti_alpha * en1kev^(1-cti_beta))) * n_transfers * trap_density_pixel * 3.65
 ;en500ev = 0.5
 ;losses500ev = (1. - exp(-cti_alpha * en500ev^(1-cti_beta))) * n_transfers * trap_density_pixel * 3.65
 
 ;losses_corrections = losses_array_ev - cti_corrections_ev_storage
 ;plot, list_energies, losses_corrections, psym=1, background=1, color=0, ytit='Losses-corrections'
 ;print, 'Mean and median of losses-corrections'
 ;print, mean(losses_corrections, /nan), median(losses_corrections)
  
; print, 'Losses at 6 keV = ' ,losses6keV
; print, 'Losses at 1 keV = ' ,losses1keV
; print, 'Losses at .5 keV = ',losses500ev
end

 
