CC = g++
CFLAGS = -Wall -c -g


INCLUDES = -I. -I$(ROOTSYS)/include -I../ -I./interface -I../CommonTools/interface

EXE_DIR = ./test/

LIB_DIR= ./lib/

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib -Wl,--no-as-needed $(ROOTGLIBS) -lHtml -lMathCore -lGenVector -lMinuit -lTMVA -lRooFit 



all: drawCherenkov drawAlignmentCorrVsRun  makeAnalysisTree alignTracking fastPlotter drawLinearity drawWaveform drawRespVsXnY drawResolution

drawHVEnergyScan: DrawTools.o
	$(CC) -Wall $(INCLUDES) -o ./bin/drawHVEnergyScan $(EXE_DIR)/drawHVEnergyScan.cpp  $(LIB_DIR)/DrawTools.o  $(ROOTFLAG) $(EXTRALIBS)


drawResolution: DrawTools.o
	$(CC) -Wall $(INCLUDES) -o ./bin/drawResolution $(EXE_DIR)/drawResolution.cpp  $(LIB_DIR)/DrawTools.o  $(ROOTFLAG) $(EXTRALIBS)



drawRespVsXnY: DrawTools.o
	$(CC) -Wall $(INCLUDES) -o ./bin/drawRespVsXnY $(EXE_DIR)/drawRespVsXnY.cpp  $(LIB_DIR)/DrawTools.o  $(ROOTFLAG) $(EXTRALIBS)

drawCherenkov:  Event.o DrawTools.o RecoTree.o
	$(CC) -Wall $(INCLUDES) -o ./bin/drawCherenkov $(EXE_DIR)/drawCherenkov.cpp $(LIB_DIR)/RecoTree.o $(LIB_DIR)/DrawTools.o  $(ROOTFLAG) $(EXTRALIBS)


makeAnalysisTree: RunHelper.o HodoCluster.o TagHelper.o EnergyCalibration.o AlignmentOfficer.o Waveform.o
	$(CC) -Wall $(INCLUDES) -o ./bin/makeAnalysisTree $(EXE_DIR)/makeAnalysisTree.cpp $(LIB_DIR)/Waveform.o $(LIB_DIR)/RunHelper.o $(LIB_DIR)/HodoCluster.o $(LIB_DIR)/TagHelper.o $(LIB_DIR)/EnergyCalibration.o $(LIB_DIR)/AlignmentOfficer.o $(ROOTFLAG) $(EXTRALIBS)

alignTracking:  HodoCluster.o AlignmentOfficer.o TagHelper.o DrawTools.o
	$(CC) -Wall $(INCLUDES) -o ./bin/alignTracking $(EXE_DIR)/alignTracking.cpp $(LIB_DIR)/HodoCluster.o $(LIB_DIR)/AlignmentOfficer.o $(LIB_DIR)/TagHelper.o $(LIB_DIR)/DrawTools.o $(ROOTFLAG) $(EXTRALIBS)

drawAlignmentCorrVsRun:	DrawTools.o
	$(CC) -Wall $(INCLUDES) -o ./bin/drawAlignmentCorrVsRun $(EXE_DIR)/drawAlignmentCorrVsRun.cpp $(LIB_DIR)/DrawTools.o $(ROOTFLAG) $(EXTRALIBS)
#analyzeCherenkov: Event.o CherenkovAnalyzer.o DrawTools.o
#	 $(CC) -Wall $(INCLUDES) -o ./bin/analyzeCherenkov $(EXE_DIR)/analyzeCherenkov.cpp  $(LIB_DIR)/Event.o $(LIB_DIR)/CherenkovAnalyzer.o  $(LIB_DIR)/DrawTools.o $(ROOTFLAG) $(EXTRALIBS)


drawWaveform: Event.o WaveformUtil.o DrawTools.o HodoCluster.o 
	 $(CC) -Wall $(INCLUDES) -o ./bin/drawWaveform $(EXE_DIR)/drawWaveform.cpp  $(LIB_DIR)/Event.o $(LIB_DIR)/WaveformUtil.o  $(LIB_DIR)/DrawTools.o  $(LIB_DIR)/HodoCluster.o $(ROOTFLAG) $(EXTRALIBS)

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

RecoTree.o: ./src/RecoTree.cc
	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/RecoTree.o -c ./src/RecoTree.cc $(ROOTFLAG) $(EXTRALIBS)


LinearityLooper.o: ./src/LinearityLooper.cc Event.o
	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/LinearityLooper.o -c ./src/LinearityLooper.cc $(ROOTFLAG) $(EXTRALIBS)

WaveformUtil.o: ./src/WaveformUtil.cc Event.o HodoCluster.o Waveform.o
	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/WaveformUtil.o $(LIB_DIR)/WaveformUtil.o -c ./src/WaveformUtil.cc $(ROOTFLAG) $(EXTRALIBS)

#CherenkovAnalyzer.o: Event.o
#	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/CherenkovAnalyzer.o -c ./src/CherenkovAnalyzer.cc $(ROOTFLAG) $(EXTRALIBS)

HodoCluster.o: ./src/HodoCluster.cc
	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/HodoCluster.o -c src/HodoCluster.cc $(ROOTFLAG) $(EXTRALIBS)

AlignmentOfficer.o: ./src/AlignmentOfficer.cc
	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/AlignmentOfficer.o -c src/AlignmentOfficer.cc $(ROOTFLAG) $(EXTRALIBS)

TagHelper.o: ./src/TagHelper.cc
	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/TagHelper.o -c src/TagHelper.cc $(ROOTFLAG) $(EXTRALIBS)

EnergyCalibration.o: ./src/EnergyCalibration.cc
	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/EnergyCalibration.o -c src/EnergyCalibration.cc $(ROOTFLAG) $(EXTRALIBS)

RunHelper.o: ../CommonTools/src/RunHelper.cc
	$(CC) $(CFLAGS) -o  $(LIB_DIR)/RunHelper.o -c ../CommonTools/src/RunHelper.cc $(ROOTFLAG) $(EXTRALIBS)

Waveform.o: ./src/Waveform.cc
	$(CC) $(CFLAGS) $(INCLUDES) -o  $(LIB_DIR)/Waveform.o -c ./src/Waveform.cc $(ROOTFLAG) $(EXTRALIBS)


clean:
	rm bin/*.o;rm lib/*.o
