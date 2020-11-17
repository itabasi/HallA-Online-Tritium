CC=gcc
CXX=g++

CFLAGS  = -O2

BINDIR = ./bin
LIBDIR = ./lib

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -lMinuit  -I /usr/local/include/eigen3  -L/usr/local/lib -lgsl -lgslcblas
ROOTGLIBS = $(shell root-config --glibs)
CXXFLAGS = -Wall -O2 $(ROOTFLAGS) 
CXXLIBS = $(ROOTLIBS)

# mp (Make Parameter)


TARGET1= xsec
OBJS1=   xsec.o Setting.o 
TARGET2= fsi
OBJS2=   fsi.o Setting.o 

.PHONY: all
all: $(TARGET1)\
     $(TARGET2)

$(LIBDIR)/%.o : %.cc
	$(CXX) $(CFLAGS) -c -o $@ $< -MMD $(CXXFLAGS)
$(TARGET1): $(patsubst %,$(LIBDIR)/%,$(OBJS1))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
$(TARGET2): $(patsubst %,$(LIBDIR)/%,$(OBJS2))
	$(CXX) $(CFLAGS) -o $(BINDIR)/$@ $^ $(CXXLIBS) $(CXXFLAGS)
.PHONY: clean
clean:
	rm -f $(LIBDIR)/*.o $(LIBDIR)/*.d core $(BINDIR)/*

