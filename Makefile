.SUFFIXES : .c .f90

CC = mpicc.mpich2
FC = mpif90.mpich2

SUBDIRS = summer_camp
CFLAGS = -g -Wall -std=c99 -pedantic
FCFLAGS = -g -Wall -pedantic
LIBS = -lm
LDFLAGS = $(LIBS)
TARGET = cpi icpi pcdf pi_montecarlo cpi2 scdf

.PHONY : clean subdirs $(SUBDIRS) all new

all : $(TARGET) subdirs

subdirs : $(SUBDIRS)

$(SUBDIRS) :
	$(MAKE) -C $@

.f90 :
	$(FC) $(FCFLAGS) -o $@ $< 

clean :
	rm -rf $(TARGET) *.o core
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir clean; done

new :
	$(MAKE) clean
	$(MAKE)
