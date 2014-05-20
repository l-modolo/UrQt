CC=g++
CFLAGS=-lz -std=c++0x -pthread -I/usr/local/include -L/usr/local/lib # -g -DEBUG -pg

all:gzstream.o Read.o Segmentation.o UrQt.cpp
	$(CC) $(CFLAGS) gzstream.o Read.o Segmentation.o UrQt.cpp -o UrQt

Segmentation.o: Segmentation.cpp
	$(CC) $(CFLAGS) -c Segmentation.cpp

Read.o: Read.cpp
	$(CC) $(CFLAGS) -c Read.cpp

gzstream.o: 
	$(CC) $(CFLAGS) -c gzstream.cpp

clean:
	rm -f *.o UrQt
