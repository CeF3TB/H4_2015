CC = g++
CFLAGS = -Wall -c -g 


SRCDIR   = src
OBJDIR   = lib
BINDIR   = bin

SOURCES  := $(wildcard $(SRCDIR)/*.cc)
BASE_INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)

INCLUDES = -I. -I$(ROOTSYS)/include -I../ -I./interface -I../CommonTools/interface 

ROOTSYS ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`

EXTRALIBS := -L$(ROOTSYS)/lib -Wl,--no-as-needed $(ROOTGLIBS) -lHtml -lMathCore -lGenVector -lMinuit -lTMVA -lRooFit `xml2-config --cflags` `xml2-config --libs`


CPP =  $(wildcard test/*.cpp)
CPPOBJ = $(patsubst test/%,%,$(CPP))
BINOBJ	=$(patsubst %.cpp,bin/%,$(CPPOBJ))
COMMONTOOLS =$(wildcard ../CommonTools/src/*.cc)
COMMONTOOLSOBJ =$(patsubst ../CommonTools/src/%.cc,lib/%.o,$(COMMONTOOLS))


all: $(BINOBJ)

$(BINOBJ) :  $(OBJECTS) $(COMMONTOOLSOBJ) 
	@echo $(BINOBJ) 	
	@echo 	@$(CC) $(INCLUDES) -o $@ $(patsubst bin/%,test/%.cpp,$@)  $(ROOTFLAG) $(EXTRALIBS) $(OBJECTS) $(COMMONTOOLSOBJ)
	@$(CC) $(INCLUDES) -o $@ $(patsubst bin/%,test/%.cpp,$@)  $(ROOTFLAG) $(EXTRALIBS) $(OBJECTS) $(COMMONTOOLSOBJ)	
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cc
	@echo $(OBJECTS)
	@echo $(CC) $(CFLAGS) $(INCLUDES) -o $@  -c $< $(ROOTFLAG) $(EXTRALIBS) 
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@  -c $< $(ROOTFLAG) $(EXTRALIBS)	
	@echo "Compiled "$<" successfully!"

$(COMMONTOOLSOBJ): 
	@echo $(COMMONTOOLSOBJ)
	@echo $(CC) $(CFLAGS)  -o $@  -c $(patsubst lib/%.o,../CommonTools/src/%.cc,$@) $(ROOTFLAG) $(EXTRALIBS)
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@  -c $(patsubst lib/%.o,../CommonTools/src/%.cc,$@) $(ROOTFLAG) $(EXTRALIBS)	
	@echo "Compiled "$<" successfully!"





clean:
	rm bin/*.o;rm lib/*.o
