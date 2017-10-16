CC = gcc
CPP = g++
CPLUSPLUS11 = -std=c++11
CFLAGS = -O3 -arch i386 -arch x86_64

all: sim

clean:
	find . -type f -name '*.o' -delete
	rm -f sim
	
sim: sim.cpp PoissonGen.o
	$(CPP) $(CFLAGS) $(CPLUSPLUS11) -o sim sim.cpp ./PoissonGen.o
	
PoissonGen.o: ../src/PoissonGen.cpp ../inc/PoissonGen.hpp
	$(CPP) $(CFLAGS) $(CPLUSPLUS11) -c -o ./PoissonGen.o ../src/PoissonGen.cpp
