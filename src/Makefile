# Makefile

OBJS = fluxRelated.o converter.o testFluxRelated.o interfaceSpeed.o minmodForVLLim.o meshClass.o configurations.o
CC = g++ 
CFLAGS = -c -Wall 
LFLAGS = 

EIGEN = -I ../dependencies/Eigen/
PYBIND = -fPIC -I /usr/include/python3.8/ $(python3 -m pybind11 --includes)


meshClass.so : meshClass.o
	$(CC) -Wall -shared -fopenmp $(PYBIND) $(EIGEN) meshClass.o -o meshClass.so

meshClass.o : meshClass.cc meshClass.h fluxRelated.h fluxRelated.cc converter.h converter.cc minmodForVLLim.h minmodForVLLim.cc interfaceSpeed.h interfaceSpeed.cc configurations.cc configurations.h
	$(CC) $(CFLAGS) -fopenmp $(PYBIND) $(EIGEN) meshClass.cc

configurations.o : configurations.cc configurations.h
	$(CC) $(CFLAGS) configurations.cc

interfaceSpeed.o : interfaceSpeed.cc interfaceSpeed.h fluxRelated.cc fluxRelated.h
	$(CC) $(CFLAGS) $(EIGEN) interfaceSpeed.cc

fluxRelated.o : fluxRelated.cc fluxRelated.h
	$(CC) $(CFLAGS) $(EIGEN) fluxRelated.cc

converter.o : converter.cc converter.h
	$(CC) $(CFLAGS) converter.cc

minmodForVLLim.o : minmodForVLLim.cc minmodForVLLim.h
	$(CC) $(CFLAGS) minmodForVLLim.cc

clean : 
	rm *.o *.so

