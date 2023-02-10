# gnu makefile for nasdaq/devel 
# builds library for nasdaq and makes executable readit
########################################################
# see http://www.gnu.org/software/make/manual/
# Here is a summary of the automatic variables.
# $@    The file name of the target.
# $%    The target member name, when the target is an archive member.
# $<    The name of the first prerequisite.
# $?    The names of all the prerequisites that are newer than the target
# $^    The names of all the prerequisites, omits duplicates
# $+                                      , retains duplicates
# $* The stem with which an implicit rule matches
########################################################
ARCH          = $(shell root-config --arch)
LOCAL        := $(shell pwd)
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
BOBJ=../bobj


CXXFLAGS     :=  $(ROOTCFLAGS) -fPIC -I$(LOCAL) 
CXXFLAGS     += -I$(LOCAL) -I/$/usr/include -I/usr/local/include -I$(BOBJ)
LIBS         := $(ROOTLIBS) L$(BOBJ) -lBaconAna  -L$(FFTW)/lib -lfftw3 $(CINTFLAGS)

ifeq ($(ARCH),macosx64)
				CXX := clang++ -arch x86_64 -ggdb 
				CXXFLAGS += -fPIC
				CPP= $(CXX) $(CXXFLAGS)
				#-install_name "@loader_path/../bobj/libBaconAna.so" -dynamiclib -o libBaconAna.so
else
				CXX := g++ -ggdb 
				CXXFLAGS += -O -Wall -fPIC 
				CPP= $(CXX) $(CXXFLAGS) 
				
endif 

#CPP+=-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64


all: sis build ana 

sis: sis.cc
	$(CPP) $< $(LIBS)  -o $@

build: build.cc
	$(CPP) $< $(LIBS)  -o $@

convertRaw: convertRaw.cc
	$(CPP) $< $(LIBS)  -o $@

ana: ana.cc anaRun.cc
	$(CPP) $< $(LIBS)  -o $@



.PHONY : clean show

clean:
	rm -f *.o 
	rm -f ana
	rm -rf ana.dSYM
show:
	@echo LD =$(LD)
	@echo sys =$(SYS)
	@echo arch =$(ARCH)
	@echo LOCAL = $(LOCAL)
	@echo BOBJ = $(BOBJ)
	@echo include path = $(INC_PATH)
	@echo CXXFLAGS = $(CXXFLAGS)
	@echo LIBS = $(LIBS)
	@echo cpp = $(CPP)
