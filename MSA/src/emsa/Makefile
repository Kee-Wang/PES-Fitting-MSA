CC := g++
CFLAG := -O2 -Wall
HEADER_PATH = ./header
INCLUDE = -I ${HEADER_PATH}


all:msa 
msa: msa.o monomial.o polynomial.o
	${CC} ${CFLAG} -o $@ $^ ${INCLUDE}; 

%:%.o
	${CC} ${CFLAG} -o $@ $^ ${INCLUDE}; 
-include .depend
depend: 
	g++ -MM *.cpp ${INCLUDE} > .depend


%.o:%.cpp
	${CC} ${CFLAG} -c $< ${INCLUDE};

# ========================================
.PHONEY:clean depend
clean:
	rm -rf *~  *.o  msa 

