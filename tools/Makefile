PBRT_TOOLS = exrassemble exravg exrtotiff tifftoexr ply2pbrt samplepat

ARCH = $(shell uname)

PBRTDIR=..
PBRTLIB=$(PBRTDIR)/core/libpbrt.a
EXRDIR=../OpenEXR
OPT=-g
CXXFLAGS = $(OPT) -I/sw/include -I$(EXRDIR)/include -I$(PBRTDIR)/core

ifeq ($(ARCH),Darwin)
LDFLAGS = -L$(EXRDIR)/lib-osx -L/sw/lib
LDFLAGS += -Wl,-all_load $(PBRTDIR)/core/libpbrt.a 
CXXFLAGS += -Wno-long-double
else
LDFLAGS = -L$(EXRDIR)/lib-linux
LDFLAGS += -Wl,--export-dynamic -Wl,-whole-archive $(PBRTLIB) -Wl,-no-whole-archive 
endif

LDFLAGS += -lIlmImf -lImath -lIex -lHalf -lz -ltiff -ldl

default: $(PBRT_TOOLS)

%: %.cpp $(PBRTLIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

ply2pbrt: ply.h ply.c ply2pbrt.c
	$(CC) $(OPT) -I. -o ply2pbrt ply2pbrt.c ply.c -lm
clean:
	rm -f $(PBRT_TOOLS) *.o *~
