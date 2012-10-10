ARCH = $(shell uname)
LEX=flex
YACC=bison -d -v -t
LEXLIB = -lfl
DLLLIB = -ldl
ifeq ($(ARCH),Darwin)
  DLLLIB =
endif
ifeq ($(ARCH),OpenBSD)
  DLLLIB =
endif

EXRINCLUDE=-I/usr/local/include/OpenEXR -I/opt/local/include/OpenEXR -I/usr/local/include/OpenEXR
EXRLIBDIR=-L/usr/local/lib -L/opt/local/lib
EXRLIBS=$(EXRLIBDIR) -Bstatic -lIex -lIlmImf -lIlmThread -lImath -lIex -lHalf -Bdynamic -lz
ifeq ($(ARCH),Linux)
  EXRLIBS += -lpthread
endif


CC=gcc
CXX=g++
LD=$(CXX) $(OPT)
DEFS=-DNDEBUG
OPT=-O2 -msse2 -mfpmath=sse
INCLUDE=-I. -Icore $(EXRINCLUDE)
WARN=-Wall
CWD=$(shell pwd)
CXXFLAGS=$(OPT) $(INCLUDE) $(WARN) $(DEFS)
CCFLAGS=$(CXXFLAGS)
LIBS=$(LEXLIB) $(DLLLIB) $(EXRLIBDIR) $(EXRLIBS) -lm 

SHARED_LDFLAGS = -shared
LRT_LDFLAGS=-rdynamic $(OPT)

ifeq ($(ARCH), Darwin)
  OS_VERSION = $(shell uname -r)
  SHARED_LDFLAGS = -flat_namespace -undefined suppress -bundle 
  LRT_LDFLAGS=$(OPT)
  INCLUDE += -I/sw/include
  #WARN += -Wno-long-double
endif

ACCELERATORS = grid kdtree
CAMERAS      = environment orthographic perspective
CORE         = api camera color dynload exrio film geometry light material mc \
               paramset parser primitive reflection sampling scene shape \
               texture timer transform transport util volume pbrtparse pbrtlex
FILM         = image
FILTERS      = box gaussian mitchell sinc triangle
INTEGRATORS  = directlighting emission irradiancecache \
               path photonmap single whitted igi debug exphotonmap
LIGHTS       = area distant goniometric infinite point projection spot infinitesample
MATERIALS    = bluepaint brushedmetal clay felt \
               glass matte mirror plastic primer \
               shinymetal skin substrate translucent uber
SAMPLERS     = bestcandidate lowdiscrepancy random stratified
SHAPES       = cone cylinder disk heightfield hyperboloid loopsubdiv nurbs \
               paraboloid sphere trianglemesh
TEXTURES     = bilerp checkerboard constant dots fbm imagemap marble mix \
               scale uv windy wrinkled
TONEMAPS     = contrast highcontrast maxwhite nonlinear
VOLUMES      = exponential homogeneous volumegrid

RENDERER     = pbrt



RENDERER_OBJS     := $(addprefix objs/, $(RENDERER:=.o) )
CORE_OBJS         := $(addprefix objs/, $(CORE:=.o) )
CORE_LIB          := core/libpbrt.a

SHAPES_DSOS       := $(addprefix bin/, $(SHAPES:=.so))
MATERIALS_DSOS    := $(addprefix bin/, $(MATERIALS:=.so))
LIGHTS_DSOS       := $(addprefix bin/, $(LIGHTS:=.so))
INTEGRATORS_DSOS  := $(addprefix bin/, $(INTEGRATORS:=.so))
VOLUMES_DSOS      := $(addprefix bin/, $(VOLUMES:=.so))
TEXTURES_DSOS     := $(addprefix bin/, $(TEXTURES:=.so))
ACCELERATORS_DSOS := $(addprefix bin/, $(ACCELERATORS:=.so))
CAMERAS_DSOS      := $(addprefix bin/, $(CAMERAS:=.so))
FILTERS_DSOS      := $(addprefix bin/, $(FILTERS:=.so))
FILM_DSOS         := $(addprefix bin/, $(FILM:=.so))
TONEMAPS_DSOS     := $(addprefix bin/, $(TONEMAPS:=.so))
SAMPLERS_DSOS     := $(addprefix bin/, $(SAMPLERS:=.so))

SHAPES_OBJS       := $(addprefix objs/, $(SHAPES:=.o))
MATERIALS_OBJS    := $(addprefix objs/, $(MATERIALS:=.o))
LIGHTS_OBJS       := $(addprefix objs/, $(LIGHTS:=.o))
INTEGRATORS_OBJS  := $(addprefix objs/, $(INTEGRATORS:=.o))
VOLUMES_OBJS      := $(addprefix objs/, $(VOLUMES:=.o))
TEXTURES_OBJS     := $(addprefix objs/, $(TEXTURES:=.o))
ACCELERATORS_OBJS := $(addprefix objs/, $(ACCELERATORS:=.o))
CAMERAS_OBJS      := $(addprefix objs/, $(CAMERAS:=.o))
FILTERS_OBJS      := $(addprefix objs/, $(FILTERS:=.o))
FILM_OBJS         := $(addprefix objs/, $(FILM:=.o))
TONEMAPS_OBJS     := $(addprefix objs/, $(TONEMAPS:=.o))
SAMPLERS_OBJS     := $(addprefix objs/, $(SAMPLERS:=.o))

RENDERER_BINARY = bin/pbrt

CORE_HEADERFILES = api.h camera.h color.h dynload.h film.h geometry.h \
                  kdtree.h light.h pbrt.h material.h mc.h mipmap.h octree.h \
                  paramset.h primitive.h reflection.h sampling.h scene.h \
                  shape.h texture.h timer.h tonemap.h transform.h transport.h \
                  volume.h 

CORE_HEADERS := $(addprefix core/, $(CORE_HEADERFILES) )

.SECONDARY: $(SHAPES_OBJS) $(MATERIALS_OBJS) $(LIGHTS_OBJS) $(INTEGRATORS_OBJS) \
            $(VOLUMES_OBJS) $(ACCELERATORS_OBJS) $(CAMERAS_OBJS) $(FILTERS_OBJS) \
            $(FILM_OBJS) $(TONEMAPS_OBJS) $(SAMPLERS_OBJS) $(TEXTURES_OBJS)

.PHONY: dirs tools exrcheck

default: dirs $(CORE_LIB) $(RENDERER_BINARY) $(INTEGRATORS_DSOS) $(VOLUMES_DSOS) $(FILM_DSOS) $(SHAPES_DSOS) $(MATERIALS_DSOS) $(LIGHTS_DSOS) $(ACCELERATORS_DSOS) $(CAMERAS_DSOS) $(SAMPLERS_DSOS) $(FILTERS_DSOS) $(TONEMAPS_DSOS) $(TEXTURES_DSOS) #tools

dirs:
	/bin/mkdir -p bin objs

tools: $(CORE_LIB)
	(cd tools && $(MAKE))

$(CORE_LIB): $(CORE_OBJS)
	@echo "Building the core rendering library (libpbrt.a)"
	@ar rcs $(CORE_LIB) $(CORE_OBJS)

bin/%.so: objs/%.o 
	@$(LD) $(SHARED_LDFLAGS) $^ -o $@

objs/%.o: renderer/%.cpp $(CORE_HEADERS)
	@echo "Building the rendering binary (pbrt)"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/%.o: core/%.cpp $(CORE_HEADERS)
	@echo "Compiling $<"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/%.o: core/%.c $(CORE_HEADERS)
	@echo "Compiling $<"
	@$(CC) $(CCFLAGS) -o $@ -c $<

objs/%.o: shapes/%.cpp $(CORE_HEADERS)
	@echo "Building Shape Plugin \"$*\""
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/%.o: integrators/%.cpp $(CORE_HEADERS)
	@echo "Building Integrator Plugin \"$*\""
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/%.o: volumes/%.cpp $(CORE_HEADERS)
	@echo "Building Volume Plugin \"$*\""
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/%.o: textures/%.cpp $(CORE_HEADERS)
	@echo "Building Texture Plugin \"$*\""
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/%.o: materials/%.cpp $(CORE_HEADERS)
	@echo "Building Material Plugin \"$*\""
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/%.o: lights/%.cpp $(CORE_HEADERS)
	@echo "Building Light Plugin \"$*\""
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/%.o: accelerators/%.cpp $(CORE_HEADERS)
	@echo "Building Accelerator Plugin \"$*\""
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/%.o: cameras/%.cpp $(CORE_HEADERS)
	@echo "Building Camera Plugin \"$*\""
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/%.o: filters/%.cpp $(CORE_HEADERS)
	@echo "Building Filter Plugin \"$*\""
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/%.o: tonemaps/%.cpp $(CORE_HEADERS)
	@echo "Building Tone Map Plugin \"$*\""
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/%.o: film/%.cpp $(CORE_HEADERS)
	@echo "Building Film Plugin \"$*\""
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

objs/%.o: samplers/%.cpp $(CORE_HEADERS)
	@echo "Building Sampler Plugin \"$*\""
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

core/pbrtlex.cpp: core/pbrtlex.l
	@echo "Lex'ing pbrtlex.l"
	@$(LEX) -o$@ core/pbrtlex.l

core/pbrtparse.h core/pbrtparse.cpp: core/pbrtparse.y
	@echo "YACC'ing pbrtparse.y"
	@$(YACC) -o $@ core/pbrtparse.y
	@if [ -e core/pbrtparse.cpp.h ]; then /bin/mv core/pbrtparse.cpp.h core/pbrtparse.h; fi
	@if [ -e core/pbrtparse.hpp ]; then /bin/mv core/pbrtparse.hpp core/pbrtparse.h; fi

$(RENDERER_BINARY): $(RENDERER_OBJS) $(CORE_LIB)
	@echo "Linking $@"
	@$(CXX) $(LRT_LDFLAGS) -o $@ $(RENDERER_OBJS) $(PBRTPRELINK) $(CORE_OBJS) $(PBRTPOSTLINK) $(LIBS)

clean:
	rm -f */*.o */*.so */*.a bin/pbrt core/pbrtlex.[ch]* core/pbrtparse.[ch]*
	(cd tools && $(MAKE) clean)

objs/exrio.o: exrcheck

exrcheck:
	@echo -n Checking for EXR installation... 
	@$(CXX) $(CXXFLAGS) -o exrcheck exrcheck.cpp $(LIBS) || \
		(cat exrinstall.txt; exit 1)
