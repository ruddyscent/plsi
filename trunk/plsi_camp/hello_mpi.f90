program  hello
include 'mpif.h'

integer ierr, rank, size

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

if (rank == 0) print *, 'first'

print *,'hello', rank, size

call MPI_FINALIZE(ierr)

end
