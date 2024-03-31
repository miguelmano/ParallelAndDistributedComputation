CC = mpicc
CFLAGS = -O2 -fopenmp

TARGET = life3d-mpi
SOURCES := $(wildcard src/*.c)
HEADERS := $(wildcard src/*.h)
OBJECTS := $(SOURCES:.c=.o)
ZIP_NAME := g27serial.zip

.PHONY: default all clean

default: $(TARGET)

all: default

$(TARGET):	$(OBJECTS)
	$(CC)	$(CFLAGS)	$(OBJECTS)	-o	$(TARGET)

%.o:	%.c
	$(CC)	$(CFLAGS)	-c	$<	-o	$@

.PRECIOUS: $(TARGET) $(OBJECTS)

clean:
	-rm -f $(OBJECTS) $(TARGET) $(ZIP_NAME)

serial:
	zip $(ZIP_NAME) $(SOURCES) $(HEADERS) Makefile