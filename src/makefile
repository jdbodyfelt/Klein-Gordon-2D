KG_OPTIONS=-DNX=256 -DNY=256 -DNZ=64
COMPILER_OPTIONS=-std=c++11 -g -fopenmp
CC=g++ $(COMPILER_OPTIONS) $(KG_OPTIONS)

INCS=
LPATH=
LIBS=

FILES=SBAB2C.cpp KG2D.cpp Tester.cpp
OBJS=$(FILES:.cpp=.o)

##################################################
Tester: clean $(OBJS)
	$(CC) $(OBJS) -o $@ $(INCS) $(LIBS) $(LPATHS) 
	make clean

%.o: %.cpp 
	$(CC) -c $< $(INCS) $(LIBS) $(LPATHS) 

clean:
	rm -rf *.o
