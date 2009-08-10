#-*- Makefile -*-
#

RHOME = d:/R/R-2.9.0

.PHONY: all libs

SOURCES=$(wildcard data/*.cpp utils/*.cpp utils/random.h model/*.cpp model/effects/*.cpp model/tables/*.cpp model/variables/*.cpp)
EXTRAS=$(foreach i,$(SOURCES),$(basename $i).o)

PKG_LIBS=-L. -lpgSn -lRMath
PKG_CPPFLAGS=-I.




all: libs SienaProfile.exe
libs: libpgSn.a

libpgSn.a: $(EXTRAS)
	 $(AR) cr libpgSn.a $(EXTRAS)

SienaProfile.exe: SienaProfile.o

clean:
	rm -f  libpgSn.a $(EXTRAS) SienaProfile.exe SienaProfile.o

CPPFLAGS = $(PKG_CPPFLAGS) -I$(RHOME)/include -DSTANDALONE
CFLAGS = -O2 -pg -DSTANDALONE
CXXFLAGS = -O2 -pg -DSTANDALONE
DLLLIBS = $(PKG_LIBS)
DLLFLAGS = -pg

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

ECHO = echo
NM = nm
DLL = g++
SED = sed

%.exe:
	$(DLL) $(DLLFLAGS) -o $@ $^ $(DLLLIBS)
