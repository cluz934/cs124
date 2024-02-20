CC=g++
CFLAGS=-O3
TARGET=randmst

all: $(TARGET)

randmst:
	$(CC) $(CFLAGS) -o $(TARGET) randmst.cpp

# Clean: remove all object files and the executable
clean:
	rm -f $(TARGET)