program reduce
include 'mpif.h'

integer ierr, myrank, nprocs, ista, iend, a(10000)
real sum, tmp

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

ista = myrank * 5000 + 1
iend = ista + 4999

do i = ista, iend
a(i) = i
enddo

sum = 0.0

do i = ista, iend
sum = sum + a(i)
enddo

call MPI_REDUCE(sum, tmp, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

sum = tmp

if (myrank == 0) then
    print *, 'sum =', sum
endif

!print *, 'rank =', myrank, ', a =', a

call MPI_FINALIZE(ierr)

end
