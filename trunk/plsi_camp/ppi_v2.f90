! ppi_v2.f90 - Fortran 90 program which calculates the value of Pi in parallel 
! written at USC parallel programming camp. This uses MPI reduce routine.
! Huioon Kim(pcandme@gist.ac.kr)

program ppi_v2
implicit none

include 'mpif.h'

integer ierr, myrank, nprocs!, idx
! integer status(MPI_STATUS_SIZE)
integer*8, parameter :: num_step = 100000
integer*8 :: i, ista, iend, stride, remain
real(kind=8) :: sum, local_pi, pi, step, x
real(kind=8) :: stime, etime

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

stride = num_step / nprocs
remain = mod(num_step, nprocs)

ista = myrank * stride + min(myrank, remain)
iend = (myrank + 1) * stride + min(myrank + 1, remain)

step = (1.0d0 / dble(num_step))
sum = 0.0d0

if (myrank .eq. 0) then
    write(*, 400)
endif

stime = MPI_WTIME()

do i = ista, iend - 1
x = (dble(i) + 0.5d0) * step
sum = sum + 4.d0 / (1.d0 + x * x)
enddo

local_pi = step * sum

call MPI_REDUCE(local_pi, pi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

etime = MPI_WTIME()

if (myrank .eq. 0) then
    write(*, 100) pi, dabs(dacos(-1.0d0) - pi)
    write(*, 300) etime - stime
    write(*, 400)
endif

100 format(' PI = ', F17.15,' (Error = ', E11.5,')')
300 format(' Elapsed Time = ', F8.4,' [sec] ')
400 format('----------------------------------------------')

call MPI_FINALIZE(ierr)

stop
end program
