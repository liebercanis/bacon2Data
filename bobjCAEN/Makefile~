# not yet ported to OSX!
########################################################
LIBNAME :=BaconAna
ARCH  = $(shell root-config --arch)
include $(ROOTSYS)/etc/Makefile.arch

SYS := $(shell uname)
LOCAL        := $(shell pwd)
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXXFLAGS     :=  $(ROOTCFLAGS) -fPIC -I$(LOCAL) 
LIBS         :=  $(ROOTGLIBS) $(CINTFLAGS)


ifeq ($(ARCH),macosx64)
				CXXFLAGS += -fPIC 
				CXX := clang++ -arch x86_64
				LD :=  clang++ -arch x86_64
				SOFLAGS := -dynamiclib -single_module

else
				LIBLINK :=
				CXXFLAGS += -fPIC 
endif 

LIBRARY := lib$(LIBNAME).so
CINTFLAGS   := -I$(ROOTSYS)/include -I/. 


#------------------------------------------------------------------------------
SKIP_OBJ = 
SKIP_DICT = 
OBJS := $(filter-out $(SKIP_OBJ),$(patsubst %.hxx,%.o,$(wildcard *.hxx)))
DICTS := $(filter-out $(SKIP_DICT), $(patsubst %.hxx,%_Dict.cxx,$(wildcard *.hxx)))
#------------------------------------------------------------------------------


all: lib 

lib: $(DICTS) $(OBJS)
	@echo "lib step $(LD) $(SOFLAGS) $(CXXFLAGS)  $(LIBS) $(OBJS) $(DICTS)  -o $(LIBRARY) ....  "
	$(LD) $(SOFLAGS) $(CXXFLAGS) $(LIBS) $(OBJS) $(DICTS) -o $(LIBRARY) 
	@echo "$@ done"


%_Dict.cxx: %.hxx %_LinkDef.h
	@echo "rootcint -f $@ -c $(CINTCFLAGS) $^"
	@rootcint -f $@ -c $(CINTFLAGS) $^
	@echo "$@ done"

	
%.o: %.cxx %.cxx
	@echo "compiling ... $(CXX) $(CXXFLAGS) "
	$(CXX) $(CXXFLAGS) -c $<
	@echo "$@ done"


.PHONY : clean show

show: 
	@echo CXX=$(CXX)
	@echo LD =$(LD)
	@echo sys =$(SYS)
	@echo arch =$(ARCH)
	@echo OBJS = $(OBJS)
	@echo DICTS = $(DICTS)
	@echo LIBS = $(LIBS)
	@echo SOFLAGS= $(SOFLAGS)
	@echo LDFLAGS= $(LDFLAGS)
	@echo CINTFLAGS= $(CINTFLAGS)
	@echo CXXFLAGS= $(CXXFLAGS)
	@echo ROOTCFLAGS= $(ROOTCFLAGS)
	@echo LIBRARY= $(LIBRARY)
	@echo LIBNAME=$(LIBNAME)	

clean:
	@echo "removing $(OBJS) *_Dict.* $(LIBRARY) *.pcm"
	@rm -f $(OBJS) *_Dict.* $(LIBRARY) *.pcm 
