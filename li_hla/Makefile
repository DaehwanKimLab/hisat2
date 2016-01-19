CXX = g++
CXXFLAGS= -W 
LINKPATH= -I./samtools-0.1.19 -L./samtools-0.1.19
LINKFLAGS = -lbam -lz -lm -lpthread 
DEBUG=
OBJECTS =

all: hla

hla: main.o 
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)

main.o: main.cpp alignments.hpp

clean:
	rm -f *.o *.gch hla
