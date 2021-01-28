# Compiler 
CC = gcc
CFLAGS = -Wall -O2 -std=gnu99

# Files
SRC = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRC))

# Name of the executable
EXECUTABLE = lmp_corr

# Compile
all: $(EXECUTABLE)

%.o: %.c
	 $(CC) $(CFLAGS) -fopenmp -c $<

# Link
$(EXECUTABLE): $(OBJS)
	 $(CC) -fopenmp $^ -o $@ -lm -lz

clean:
	 rm *.o $(EXECUTABLE)
