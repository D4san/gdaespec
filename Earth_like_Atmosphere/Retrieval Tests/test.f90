program test
  open(unit=55, file='contam_10T_0.26spot-0.70fac-.txt', status='unknown')
  write(55, *) 'Hello'
  close(55)
end program test
