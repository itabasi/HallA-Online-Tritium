CC=gcc
CXX=g++

CFLAGS  = -O2

BINDIR = ./bin
LIBDIR = ./lib

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lMinuit
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

# mp (Make Parameter)
#TARGET1= hrs_tof_twcorr
#OBJS1=   hrs_tof_twcorr.o 

TARGET1= coincalib
OBJS1=   coincalib.o Setting.o Tree.o ParamMan.o
TARGET2= t0corr
OBJS2=   t0corr.o Setting.o Tree.o ParamMan.o
TARGET3= t0calib
OBJS3=   t0calib.o Setting.o Tree.o ParamMan.o
TARGET4= t0calib_hist
OBJS4=   t0calib_hist.o Setting.o Tree.o ParamMan.o

.PHONY: all
all: $(TARGET1)\
     $(TARGET2)\
     $(TARGET3)\
     $(TARGET4)
$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< -MMD $(CXXFLAGS)

$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET2): $(patsubst %,$(LIBDIR)/%,$(OBJS2))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET3): $(patsubst %,$(LIBDIR)/%,$(OBJS3))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET4): $(patsubst %,$(LIBDIR)/%,$(OBJS4))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o $(LIBDIR)/*.d core $(BINDIR)/*

