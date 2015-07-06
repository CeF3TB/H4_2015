CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include -I../ -I./interface -I../CommonTools/interface

EXE_DIR = ./test/

LIB_DIR= ./lib/

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib -Wl,--no-as-needed $(ROOTGLIBS) -lHtml -lMathCore -lGenVector -lMinuit -lTMVA -lRooFit 



all: fastPlotter


fastPlotter: fastPlotter fastPlotter.o DrawTools.o
	 $(CC) -Wall $(INCLUDES) -o ./bin/fastPlotter $(EXE_DIR)/fastPlotter.cpp $(LIB_DIR)/fastPlotter.o  $(LIB_DIR)/DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

DrawTools.o: ../CommonTools/src/DrawTools.cc
	$(CC) $(CFLAGS) -o  $(LIB_DIR)/DrawTools.o -c ../CommonTools/src/DrawTools.cc $(ROOTFLAG) $(EXTRALIBS)

fastPlotter.o: ./src/fastPlotter.cc
	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/fastPlotter.o -c ./src/fastPlotter.cc $(ROOTFLAG) $(EXTRALIBS)




clean:
	rm *.o
