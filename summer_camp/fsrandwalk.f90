program serial_randwalk
implicit none

integer :: i, j, min_p, max_p
integer, parameter :: N_P = 100000, N_W = 100000
integer, dimension(N_P) :: part
integer, dimension(:), allocatable :: pos
real :: rand, temp

part = 0

call srand(1)

do i = 1, N_P
do j = 1, N_W
temp = rand(0)
if (temp <= 0.5) then
    part(i) = part(i) + 1
else
    part(i) = part(i) - 1
endif
enddo
enddo

min_p = minval(part)
max_p = maxval(part)
allocate(pos(min_p:max_p))

do i = 1, N_P
pos(part(i)) = pos(part(i)) + 1
enddo

do i = min_p, max_p, 10
print *, i, pos(i)
enddo

deallocate(pos)

end program serial_randwalk
