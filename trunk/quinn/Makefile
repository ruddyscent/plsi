.SUFFIXES : .c .f90

CC = mpicc.mpich2
FC = mpif90.mpich2

CFLAGS = -g -Wall -std=c99 -pedantic
FCFLAGS = -g -Wall
LIBS = -lm
LDFLAGS = $(LIBS)

TARGET = fig8-8 fig8-14 gen_matrix gen_vector

.PHONY : clean cleandiff cleanobj cleanbackup all new

all : $(TARGET)

.f90:
	$(FC) $(FCFLAGS) -o $@ $< 

clean : cleandiff cleanobj cleanbackup
	rm -f $(TARGET)

leandiff :
	rm -f *.diff

cleanobj :
	rm -f *.o

cleanbackup :
	rm -f *~

new :
	$(MAKE) clean
	$(MAKE)
