.SUFFIXES : .c .o

CC = mpicc

INCLUDE = 
LIBS = -lm
CFLAGS = -g -Wall -std=c99 $(INCLUDE)

TARGET = cpi icpi cdf

all : $(TARGET)

clean :
	rm -rf $(OBJS) $(TARGET) core

new :
	$(MAKE) clean
	$(MAKE)
