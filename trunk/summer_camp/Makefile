.SUFFIXES : .c .f90

CC = mpicc.mpich2
FC = mpif90.mpich2

CFLAGS = -g -Wall -std=c99 -pedantic
FCFLAGS = -g -Wall
LIBS = -lm
LDFLAGS = $(LIBS)

TARGET = fppi_v1 fppi_v2 send isend bcast gather gatherv allgather reduce scatter fgather fprandwalk fspi hello fbcast fisend freduce fsrandwalk hello_mpi

.PHONY : clean cleandiff cleanobj cleanbackup all new

all : $(TARGET)

.f90:
	$(FC) $(FCFLAGS) -o $@ $< 

clean : cleandiff cleanobj cleanbackup
	rm -f $(TARGET)

cleandiff :
	rm -f *.diff

cleanobj :
	rm -f *.o

cleanbackup :
	rm -f *~

new :
	$(MAKE) clean
	$(MAKE)
