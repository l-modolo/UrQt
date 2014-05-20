all:gzstream.cpp Read.cpp polyNtrimmer.cpp
	g++ gzstream.cpp Read.cpp polyNtrimmer.cpp -o polyNtrimmer -lz -std=c++0x -pthread -I/usr/local/include -L/usr/local/lib# -g -DEBUG -pg -Werror

clean:
	rm -f *.o polyAtrimmer
