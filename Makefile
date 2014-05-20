CC=g++
CFLAGS=-lz -std=c++0x -pthread -I/usr/local/include -L/usr/local/lib -g -DEBUG -pg -Werror

all:gzstream.o Read.o Segmentation.o polyNtrimmer.cpp
	$(CC) $(CFLAGS) gzstream.o Read.o Segmentation.o polyNtrimmer.cpp -o polyNtrimmer

Segmentation.o: Segmentation.cpp
	$(CC) $(CFLAGS) -c Segmentation.cpp

Read.o: Read.cpp
	$(CC) $(CFLAGS) -c Read.cpp

gzstream.o: 
	rm -f *.o polyAtrimmer
	$(CC) $(CFLAGS) -c gzstream.cpp

clean:
	rm -f *.o polyAtrimmer
