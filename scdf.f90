! scdf.f90 - Fortran 90 serial code which calculates the value of CDF of the standard normal distribution in serial order.

program scdf
implicit none

integer, parameter :: NUM_STEP = 100000
integer :: i
real(8) :: p = 0.5d0, x = 0.0d0, z = 1.0d0

do i = 0, NUM_STEP
p = p + 1.0d0 / sqrt(2.0 * acos(-1.0d0)) * exp(-x ** 2 / 2.0d0) * z / real(NUM_STEP, 8)
x = x + z / real(NUM_STEP, 8)
enddo

print *, 'The probability is ', p, ' Z = ', z

end program scdf
