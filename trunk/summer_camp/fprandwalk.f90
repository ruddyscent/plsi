program fprandwalk
implicit none

include 'mpif.h'

integer :: i, j, d, min_p, max_p, ierr, myrank, nprocs
integer :: ista, iend
integer, parameter :: N_P = 100000, N_W = 100000
integer, dimension(N_P) :: part, tmp
integer, dimension(:), allocatable :: pos
real :: rand, temp

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

ista = myrank * (N_P / nprocs) + 1
iend = ista + (N_P / nprocs) - 1

part = 0
tmp = 0

call srand(0.5 + myrank * 0.1)

do i = ista, iend
do j = 1, N_W
temp = rand(0)
if (temp <= 0.5) then
    part(i) = part(i) + 1
else
    part(i) = part(i) - 1
endif
enddo
enddo

call MPI_REDUCE(part, tmp, N_P, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

part = tmp

if (myrank .eq. 0) then
    min_p = minval(part)
    max_p = maxval(part)
    allocate(pos(min_p:max_p))
    pos = 0

    do i = 1, N_P
    pos(part(i)) = pos(part(i)) + 1
    enddo

    do i = min_p, max_p, 10
    print *, i, pos(i)
    enddo

    deallocate(pos)
endif

call MPI_FINALIZE(ierr)

end program parallel_randwalk
