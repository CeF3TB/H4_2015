CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include -I../ -I./interface -I../CommonTools/interface

EXE_DIR = ./test/

LIB_DIR= ./lib/

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib -Wl,--no-as-needed $(ROOTGLIBS) -lHtml -lMathCore -lGenVector -lMinuit -lTMVA -lRooFit 



all: fastPlotter drawLinearity drawWaveform

drawWaveform: Event.o WaveFormUtil.o DrawTools.o
	 $(CC) -Wall $(INCLUDES) -o ./bin/drawWaveform $(EXE_DIR)/drawWaveform.cpp  $(LIB_DIR)/Event.o $(LIB_DIR)/WaveFormUtil.o  $(LIB_DIR)/DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawLinearity: Event.o LinearityLooper.o DrawTools.o
	 $(CC) -Wall $(INCLUDES) -o ./bin/drawLinearity $(EXE_DIR)/drawLinearity.cpp $(LIB_DIR)/Event.o $(LIB_DIR)/LinearityLooper.o  $(LIB_DIR)/DrawTools.o $(ROOTFLAG) $(EXTRALIBS)	

fastPlotter: fastPlotter.o DrawTools.o
	 $(CC) -Wall $(INCLUDES) -o ./bin/fastPlotter $(EXE_DIR)/fastPlotter.cpp $(LIB_DIR)/fastPlotter.o  $(LIB_DIR)/DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

DrawTools.o: ../CommonTools/src/DrawTools.cc
	$(CC) $(CFLAGS) -o  $(LIB_DIR)/DrawTools.o -c ../CommonTools/src/DrawTools.cc $(ROOTFLAG) $(EXTRALIBS)

fastPlotter.o: ./src/fastPlotter.cc
	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/fastPlotter.o -c ./src/fastPlotter.cc $(ROOTFLAG) $(EXTRALIBS)

Event.o: ./src/Event.cc
	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/Event.o -c ./src/Event.cc $(ROOTFLAG) $(EXTRALIBS)

LinearityLooper.o: ./src/LinearityLooper.cc Event.o
	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/LinearityLooper.o -c ./src/LinearityLooper.cc $(ROOTFLAG) $(EXTRALIBS)

WaveFormUtil.o: ./src/WaveFormUtil.cc Event.o
	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/WaveFormUtil.o -c ./src/WaveFormUtil.cc $(ROOTFLAG) $(EXTRALIBS)


clean:
	rm bin/*.o;rm lib/*.o
