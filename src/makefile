KG_OPTIONS=-DNX=512 -DNY=512 -DNZ=128
COMPILER_OPTIONS=-std=c++11 -g 
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
