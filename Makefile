CC=gcc
CFLAGS=-std=c++0x -pthread -I/usr/local/include -L/usr/lib -lz
CFLAGSTATIC=-std=c++0x -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -I/usr/local/include -L/usr/lib -lz -static-libgcc -static-libstdc++ -static

all:gzstream.o Read.o Segmentation.o UrQt.cpp
	$(CC) $(CFLAGS) gzstream.o Read.o Segmentation.o UrQt.cpp -o UrQt

Segmentation.o: Segmentation.cpp
	$(CC) $(CFLAGS) -c Segmentation.cpp

Read.o: Read.cpp
	$(CC) $(CFLAGS) -c Read.cpp

gzstream.o: 
	$(CC) $(CFLAGS) -c gzstream.cpp

static:gzstream_static.o Read_static.o Segmentation_static.o UrQt.cpp
	$(CC) $(CFLAGSTATIC) gzstream.o Read.o Segmentation.o UrQt.cpp /usr/lib/libz.a -o UrQt

Segmentation_static.o: Segmentation.cpp
	$(CC) $(CFLAGSTATIC) -c Segmentation.cpp

Read_static.o: Read.cpp
	$(CC) $(CFLAGSTATIC) -c Read.cpp

gzstream_static.o: 
	$(CC) $(CFLAGSTATIC) -c gzstream.cpp

clean:
	rm -f *.o UrQt