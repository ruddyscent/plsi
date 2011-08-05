.SUFFIXES : .c .o

CC = mpicc

INCLUDE = 
CFLAGS = -g -Wall -std=c99 $(INCLUDE)
LIBS = -lm
LDFLAGS = $(LIBS)

TARGET = cpi icpi cdf

all : $(TARGET)

clean :
	rm -rf $(OBJS) $(TARGET) core

new :
	$(MAKE) clean
	$(MAKE)
