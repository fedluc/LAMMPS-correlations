# Compiler 
CC = gcc
CFLAGS = -Wall -O2 -std=gnu99

# Files
SRC = $(wildcard *.c)
OBJS = $(patsubst %.c, %.o, $(SRC))

# Name of the executable
EXECUTABLE = isf

# Compile
all: $(EXECUTABLE)

%.o: %.c
	 $(CC) $(CFLAGS) $(INCLUDE) -c $<

# Link
$(EXECUTABLE): $(OBJS)
	 $(CC) $(LIB) $^ -o $@ -lz

clean:
	 rm *.o isf
