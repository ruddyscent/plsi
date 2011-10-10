.SUFFXES : .c .f90

CC = mpicc
FC = mpif90

CFLAGS = -g -Wall -O0 -std=c99
LIBS = -lm
LDFLAGS = $(LIBS)
FCFLAGS = -g -Wall -O0
TARGET = problem2_1 problem2_2 problem2_3 problem2_1_mpi problem2_2_mpi problem2_3_mpi

.PHONY : clean cleandiff cleanobj cleanbackup all new

all : $(TARGET)

.f90 :
	$(FC) $(FCFLAGS) -o $@ $<

clean : cleandiff cleanobj cleanbackup
	rm -f $(TARGET)

cleandiff :
	rm -f *.diff

cleanobj :
	rm -f *.o

cleanbackup :
	rm -f *~

cleanswap :
	rm -f \#*\#

new :
	$(MAKE) clean
	$(MAKE)
