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
ROOTINCS     := $(shell root-config --incdir)
XCODEDIR =/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1
XCODEDIR2=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include/c++/v1/
ROOT_INCLUDE_PATH += $(XCODEDIR) $(XCODEDIR2)

CXXFLAGS     :=  $(ROOTCFLAGS) -fPIC -I$(LOCAL) 
LIBS         :=  $(ROOTLIBS) $(CINTFLAGS)


ifeq ($(ARCH),macosx64)
				CXX := clang++ -arch x86_64 $(CXXFLAGS)
				LD :=  clang++ -arch x86_64
				SOFLAGS := -dynamiclib -single_module

else
				LIBLINK :=
				CXXFLAGS += -fPIC 
endif 

LIBRARY := lib$(LIBNAME).so
CINTFLAGS   := -I$(ROOTINCS) 
#-I/usr/local/Cellar/llvm/16.0.1/include/c++/v1 -I/usr/local/Cellar/llvm/16.0.1/./lib/clang/16/include -I/usr/local/Cellar/llvm/16.0.1/./Toolchains/LLVM16.0.1.xctoolchain/usr/include -I/usr/local/Cellar/gcc/12.2.0 
#-I$(XCODEDIR) -I$(XCODEDIR2) -I/. 


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
	@echo "rootcint -f  $@ -c $(CINTFLAGS) $^"
	@rootcint -f   $@ -c $(CINTFLAGS) $^
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
	@echo ROOTINCS= $(ROOTINCS)
	@echo LIBRARY= $(LIBRARY)
	@echo LIBNAME=$(LIBNAME)	

clean:
	@echo "removing $(OBJS) *_Dict.* $(LIBRARY) *.pcm"
	@rm -f $(OBJS) *_Dict.* $(LIBRARY) *.pcm 
