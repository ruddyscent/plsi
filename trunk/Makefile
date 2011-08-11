.SUFFIXES : .c .f90

CC = mpicc
FC = mpif90

CFLAGS = -g -Wall -std=c99 -pedantic
FCFLAGS = -g -Wall 
LIBS = -lm
LDFLAGS = $(LIBS)

TARGET = cpi icpi cdf pi_montecarlo cpi2 scdf plsi_camp/ppi_v1 plsi_camp/ppi_v2 plsi_camp/send plsi_camp/isend plsi_camp/bcast plsi_camp/gather plsi_camp/gatherv plsi_camp/allgather plsi_camp/reduce plsi_camp/scatter

all : $(TARGET)

.f90:
	$(FC) $(FCFLAGS) -o $@ $< 

clean :
	rm -rf $(OBJS) $(TARGET) *.o core

clean_all :
	rm -rf $(OBJS) $(TARGET) *.o core *~

new :
	$(MAKE) clean
	$(MAKE)
