program parallel_pi_v2
implicit none

include 'mpif.h'

integer ierr, idx, myrank, nprocs
integer status(MPI_STATUS_SIZE)
integer*8, parameter :: num_step = 5000000000
integer*8 :: i, ista, iend
real(kind=8) :: sum, tmp, step, pi, x

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

ista = myrank * (num_step / nprocs) + 1
iend = ista + (num_step / nprocs) - 1

step = (1.0d0 / dble(num_step))
sum = 0.0d0

if (myrank .eq. 0) then
    write(*, 400)
endif

do i = ista, iend 
x = (dble(i) - 0.5d0) * step
sum = sum + 4.d0 / (1.d0 + x * x)
enddo

pi = step * sum

call MPI_REDUCE(pi, tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

pi = tmp

if (myrank .eq. 0) then
    write(*, 100) pi, dabs(dacos(-1.0d0) - pi)
    write(*, 400)
endif

100 format(' PI = ', F17.15,' (Error = ', E11.5,')')
400 format('----------------------------------------------')
stop

call MPI_FINALIZE(ierr)

end program
