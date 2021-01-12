#Makefile for testing cause I am too pissed off 

OBJS = fluxRelated.o converter.o testFluxRelated.o interfaceSpeed.o minmodForVLLim.o
CC = g++ -std=c++14
CFLAGS = -c -Wall 
LFLAGS = 
INCLUDE = -I ../dependencies/Eigen/


all : $(OBJS) 
	$(CC) $(LFLAGS) $(OBJS) -o exec 

testFluxRelated.o : testFluxRelated.cc fluxRelated.cc converter.cc fluxRelated.h converter.h
	$(CC) $(CFLAGS) $(INCLUDE) testFluxRelated.cc

interfaceSpeed.o : interfaceSpeed.cc interfaceSpeed.h fluxRelated.cc fluxRelated.h
	$(CC) $(CFLAGS) $(INCLUDE) interfaceSpeed.cc

fluxRelated.o : fluxRelated.cc fluxRelated.h
	$(CC) $(CFLAGS) $(INCLUDE) fluxRelated.cc

converter.o : converter.cc converter.h
	$(CC) $(CFLAGS) converter.cc

minmodForVLLim.o : minmodForVLLim.cc minmodForVLLim.h
	$(CC) $(CFLAGS) minmodForVLLim.cc
clean : 
	\rm *.o exec

