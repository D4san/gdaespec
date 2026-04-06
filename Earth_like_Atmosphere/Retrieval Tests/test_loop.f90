program test
  integer i
  do i = 1, 100000
    open(unit=55, file='contam_10T_0.26spot-0.70fac-.txt', status='replace')
    write(55, *) 'Hello', i
    close(55)
  end do
end program test
