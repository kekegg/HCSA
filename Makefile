#############################################################################
# Makefile for building: efba
# Generated by qmake (2.01a) (Qt 4.4.2) on: Tue Sep 23 14:17:57 2008
# Project:  efba.pro
# Template: app
# Command: /usr/bin/qmake -unix -o Makefile efba.pro
#############################################################################

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = -DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED
CFLAGS        = -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)
CXXFLAGS      = -g -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)
INCPATH       = -I/usr/share/qt4/mkspecs/linux-g++ -I. -I. -I. -I.
LINK          = g++
LFLAGS        = -Wl
LIBS          = $(SUBLIBS)  -L/usr/lib -lpthread -lglpk -lgsl -lgslcblas
AR            = ar cqs
RANLIB        = 
QMAKE         = /usr/bin/qmake
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = $(COPY)
COPY_DIR      = $(COPY) -r
INSTALL_FILE  = install -m 644 -p
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = install -m 755 -p
DEL_FILE      = rm -f
SYMLINK       = ln -sf
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = func.cpp \
		hcsa.cpp \
		main.cpp 
OBJECTS       = func.o \
		hcsa.o \
		main.o

QMAKE_TARGET  = hcsa
DESTDIR       = 
TARGET        = hcsa

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: $(TARGET)

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

qmake:  FORCE
	@$(QMAKE) -unix -o Makefile efba.pro


clean:compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Compile

func.o: func.cpp func.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o func.o func.cpp

hcsa.o: hcsa.cpp func.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o hcsa.o hcsa.cpp

main.o: main.cpp func.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:
