KG_OPTIONS=-DNX=256 -DNY=256 -DNZ=2
COMPILER_OPTIONS=-std=c++11 -g 
CC=g++ $(COMPILER_OPTIONS) $(KG_OPTIONS)

AF_PATH=/usr/local
INCS=#-I$(AF_PATH)/include
LPATH=#-L$(AF_PATH)/lib
LIBS=#-lafcpu

FILES=SBAB2C.cpp KG2D.cpp Tester.cpp
OBJS=$(FILES:.cpp=.o)

##################################################
Tester: clean $(OBJS)
	$(CC) $(OBJS) -o $@ $(INCS) $(LIBS) $(LPATHS) 

%.o: %.cpp 
	$(CC) -c $< $(INCS) $(LIBS) $(LPATHS) 

clean:
	rm -rf *.o
