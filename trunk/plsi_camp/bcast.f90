program bcast
include 'mpif.h'

integer i, ierr, rank, size, imsg(4)

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

if (rank == 0) then
    do i = 1, 4
    imsg(i) = i
    enddo
else
    do i = 1, 4
    imsg(i) = 0
    enddo
endif

print *, 'rank:', rank, ', before:', imsg

call MPI_BCAST(imsg, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

print *, 'rank:', rank, ', after: ', imsg

call MPI_FINALIZE(ierr)

end
