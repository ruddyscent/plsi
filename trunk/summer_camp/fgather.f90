program fgather
include 'mpif.h'

integer ierr, myrank, nprocs, irecv(3)

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

isend = myrank + 1

call MPI_GATHER(isend, 1, MPI_INTEGER, irecv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

if (myrank == 0) then
    print *, 'irecv =', irecv
endif

call MPI_FINALIZE(ierr)

end
