program fisend
include 'mpif.h'

integer err, rank, size, count
real data(100), value(200)
integer status(MPI_STATUS_SIZE)

call MPI_INIT(err)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, err)

if (rank .eq. 0) then
    data = 3.0
    call MPI_SEND(data, 100, MPI_REAL, 1, 55, MPI_COMM_WORLD, err)
elseif (rank .eq. 1) then
    call MPI_RECV(value, 200, MPI_REAL, MPI_ANY_SOURCE, 55, MPI_COMM_WORLD, status, err)
    print *, "P:", rank, " got data from processor ", status(MPI_SOURCE)
    call MPI_GET_COUNT(status, MPI_REAL, count, err)
    print *, "P:", rank, " got ", count, " elements"
    print *, "P:", rank, " value(5)=", value(5)
endif
call MPI_FINALIZE(err)
end
