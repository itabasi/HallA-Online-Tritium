#Makefile

OBJS = main.o CalcKinematics_eek.o
TARGET = KINEMA_EEK

#ROOTGLIBS = -lGpad -lGui -lGX11 -lX3d
#ROOTLIBS = `root-config --libs`
#ROOTLIBS  = -L/cern/root_v5.24.00.x86_64_fc8_gcc4.1.2/lib/ -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -pthread -lm -ldl -rdynamic
#ROOTLIBS = -L/home/gogami/root/root-5.34.36-install/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic
#ROOTINCLUDE = -pthread -m64 -I/home/gogami/root/root-5.34.36-install/include/root
ROOTLIBS = -L/usr/local/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
ROOTINCLUDE = -pthread -std=c++11 -m64 -I/usr/local/include
#ROOTINCLUDE = `root-config --cflags`
GLIBS = $(ROOTLIBS) -lm 

CC = g++
#CC = icc

CCFLAGS = $(ROOTINCLUDE)
#CCFLAGS = -O  -D__osf__ -D__alpha -I$(ROOTINCLUDE)
#LD = g++
#LDFLAGS = -g

all: ${TARGET}

${TARGET}: ${OBJS}
	${CC} ${CCFLAGS} -o $@ ${OBJS} ${GLIBS}

.cc.o:
	${CC} ${CCFLAGS} -c $<

clean:
	rm -rf *.o core ${TARGET}

