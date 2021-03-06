.SUFFIXES : .c .f90

CC = mpicc
FC = mpif90

SUBDIRS = summer_camp quinn main
CFLAGS = -g -Wall -std=c99 -pedantic
FCFLAGS = -g -Wall
LIBS = -lm
LDFLAGS = $(LIBS)
TARGET = cpi icpi pcdf pi_montecarlo cpi2 scdf pio_ex1 pio_ex2 pio_ex3 pio_ex4 mvmult

.PHONY : clean cleandiff cleanobj cleanbackup subdirs $(SUBDIRS) all new

all : $(TARGET) subdirs

subdirs : $(SUBDIRS)

$(SUBDIRS) :
	$(MAKE) -C $@

.f90 :
	$(FC) $(FCFLAGS) -o $@ $< 

clean : cleandiff cleanobj cleanbackup
	rm -f $(TARGET)
	$(MAKE) -C $(SUBDIRS) clean

cleandiff :
	rm -f *.diff

cleanobj :
	rm -f *.o

cleanbackup :
	rm -f *~

new :
	$(MAKE) clean
	$(MAKE)
