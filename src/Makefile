# Makefile

OBJS = fluxRelated.o converter.o interfaceSpeed.o minmodForVLLim.o meshClass.o configurations.o
CC = g++ 
CFLAGS = -c -Wall 
LFLAGS = 

EIGEN = -I ../dependencies/Eigen/
PYBIND = -fPIC -I /usr/include/python3.8/ $(python3 -m pybind11 --includes)


meshClass.so : meshClass.o configurations.o interfaceSpeed.o fluxRelated.o converter.o minmodForVLLim.o
	$(CC) -O3 -Wall -shared -fopenmp $(PYBIND) $(EIGEN) $(OBJS) -o meshClass.so

meshClass.o : meshClass.cc meshClass.h fluxRelated.h converter.h minmodForVLLim.h interfaceSpeed.h configurations.h
	$(CC) $(CFLAGS) -fopenmp $(PYBIND) $(EIGEN) meshClass.cc

configurations.o : configurations.cc configurations.h converter.h
	$(CC) $(CFLAGS) configurations.cc

interfaceSpeed.o : interfaceSpeed.cc interfaceSpeed.h fluxRelated.h
	$(CC) $(CFLAGS) $(EIGEN) interfaceSpeed.cc

fluxRelated.o : fluxRelated.cc fluxRelated.h
	$(CC) $(CFLAGS) $(EIGEN) fluxRelated.cc

converter.o : converter.cc converter.h
	$(CC) $(CFLAGS) converter.cc

minmodForVLLim.o : minmodForVLLim.cc minmodForVLLim.h
	$(CC) $(CFLAGS) minmodForVLLim.cc

clean : 
	rm *.o *.so

cleanData : 
	rm *.csv
