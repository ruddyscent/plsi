.SUFFIXES : .c .f90

CC = mpicc
FC = mpif90

CFLAGS = -g -Wall -std=c99
FFLAGS = -g -Wall
LIBS = -lm
LDFLAGS = $(LIBS)

TARGET = cpi icpi cdf pi_montecarlo cpi2 ppi_v1 ppi_v2

all : $(TARGET)

.f90:
	$(FC) $(FFLAGS) -o $@ $< 

clean :
	rm -rf $(OBJS) $(TARGET) *.o core

clean_all :
	rm -rf $(OBJS) $(TARGET) *.o core *~

new :
	$(MAKE) clean
	$(MAKE)
