
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

// api.cpp*
#include "api.h"
#include "paramset.h"
#include "color.h"
#include "scene.h"
#include "film.h"
#include "dynload.h"
#include "volume.h"
#include <map>
using std::map;
#if (_MSC_VER >= 1400) // NOBOOK
#include <stdio.h>     // NOBOOK
#define snprintf _snprintf // NOBOOK
#endif // NOBOOK
// API Local Classes
struct RenderOptions {
	// RenderOptions Public Methods
	RenderOptions();
	Scene *MakeScene() const;
	// RenderOptions Public Data
	string FilterName;
	ParamSet FilterParams;
	string FilmName;
	ParamSet FilmParams;
	string SamplerName;
	ParamSet SamplerParams;
	string AcceleratorName;
	ParamSet AcceleratorParams;
	string SurfIntegratorName, VolIntegratorName;
	ParamSet SurfIntegratorParams, VolIntegratorParams;
	string CameraName;
	ParamSet CameraParams;
	Transform WorldToCamera;
	bool gotSearchPath;
	mutable vector<Light *> lights;
	mutable vector<Reference<Primitive> > primitives;
	mutable vector<VolumeRegion *> volumeRegions;
	map<string, vector<Reference<Primitive> > > instances;
	vector<Reference<Primitive> > *currentInstance;
};
RenderOptions::RenderOptions() {
	// RenderOptions Constructor Implementation
	FilterName = "mitchell";
	FilmName = "image";
	SamplerName = "bestcandidate";
	AcceleratorName = "kdtree";
	SurfIntegratorName = "directlighting";
	VolIntegratorName = "emission";
	CameraName = "perspective";
	char *searchEnv = getenv("PBRT_SEARCHPATH");
	if (searchEnv == NULL) {
		Warning("PBRT_SEARCHPATH not set in your environment.\n"
			  "pbrt won't be able to find plugins if "
			  "no SearchPath in input file.\n"
			  "PBRT_SEARCHPATH should be a "
			  "\"%s\"-separated list of directories.\n",
			  PBRT_PATH_SEP);
		gotSearchPath = false;
	}
	else {
		UpdatePluginPath(searchEnv);
		gotSearchPath = true;
	}
	currentInstance = NULL;
}
struct GraphicsState {
	// Graphics State Methods
	GraphicsState();
	// Graphics State
	map<string, Reference<Texture<float> > >
		floatTextures;
	map<string, Reference<Texture<Spectrum> > >
		spectrumTextures;
	ParamSet materialParams;
	string material;
	ParamSet areaLightParams;
	string areaLight;
	bool reverseOrientation;
};
GraphicsState::GraphicsState() {
	// GraphicsState Constructor Implementation
	material = "matte";
	reverseOrientation = false;
}
// API Static Data
#define STATE_UNINITIALIZED  0
#define STATE_OPTIONS_BLOCK  1
#define STATE_WORLD_BLOCK    2
static int currentApiState = STATE_UNINITIALIZED;
static Transform curTransform;
static map<string, Transform> namedCoordinateSystems;
static RenderOptions *renderOptions = NULL;
static GraphicsState graphicsState;
static vector<GraphicsState> pushedGraphicsStates;
static vector<Transform> pushedTransforms;
// API Macros
#define VERIFY_INITIALIZED(func) \
if (currentApiState == STATE_UNINITIALIZED) { \
	Error("pbrtInit() must be before calling \"%s()\". " \
		"Ignoring.", func); \
	return; \
} else /* swallow trailing semicolon */
#define VERIFY_OPTIONS(func) \
VERIFY_INITIALIZED(func); \
if (currentApiState == STATE_WORLD_BLOCK) { \
	Error("Options cannot be set inside world block; " \
		"\"%s\" not allowed.  Ignoring.", func); \
	return; \
} else /* swallow trailing semicolon */
#define VERIFY_WORLD(func) \
VERIFY_INITIALIZED(func); \
if (currentApiState == STATE_OPTIONS_BLOCK) { \
	Error("Scene description must be inside world block; " \
		"\"%s\" not allowed. Ignoring.", func); \
	return; \
} else /* swallow trailing semicolon */
// API Function Definitions
COREDLL void pbrtInit() {
	// System-wide initialization
	// Make sure floating point unit's rounding stuff is set
	// as is expected by the fast FP-conversion routines.  In particular,
	// we want double precision on Linux, not extended precision!
	#ifdef FAST_INT
	#if defined(__linux__) && defined(__i386__)
	int cword = _FPU_MASK_DM | _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_PM |
		_FPU_MASK_UM | _FPU_MASK_IM | _FPU_DOUBLE | _FPU_RC_NEAREST;
	_FPU_SETCW(cword);
	#endif
	#if defined(WIN32)
	_control87(_PC_53, MCW_PC);
	#endif
	#endif // FAST_INT
	// API Initialization
	if (currentApiState != STATE_UNINITIALIZED)
		Error("pbrtInit() has already been called.");
	currentApiState = STATE_OPTIONS_BLOCK;
	renderOptions = new RenderOptions;
	graphicsState = GraphicsState();
}
COREDLL void pbrtCleanup() {
	StatsCleanup();
	// API Cleanup
	if (currentApiState == STATE_UNINITIALIZED)
		Error("pbrtCleanup() called without pbrtInit().");
	else if (currentApiState == STATE_WORLD_BLOCK)
		Error("pbrtCleanup() called while inside world block.");
	currentApiState = STATE_UNINITIALIZED;
	delete renderOptions;
	renderOptions = NULL;
}
COREDLL void pbrtIdentity() {
	VERIFY_INITIALIZED("Identity");
	curTransform = Transform();
}
COREDLL void pbrtTranslate(float dx, float dy, float dz) {
	VERIFY_INITIALIZED("Translate");
	curTransform =
		curTransform * Translate(Vector(dx, dy, dz));
}
COREDLL void pbrtTransform(float tr[16]) {
	VERIFY_INITIALIZED("Transform");
	curTransform = Transform(new Matrix4x4(
		tr[0], tr[4], tr[8], tr[12],
		tr[1], tr[5], tr[9], tr[13],
		tr[2], tr[6], tr[10], tr[14],
		tr[3], tr[7], tr[11], tr[15]));
}
COREDLL void pbrtConcatTransform(float tr[16]) {
	VERIFY_INITIALIZED("ConcatTransform");
	curTransform = curTransform * Transform(
		new Matrix4x4(tr[0], tr[4], tr[8], tr[12],
				 tr[1], tr[5], tr[9], tr[13],
				 tr[2], tr[6], tr[10], tr[14],
				 tr[3], tr[7], tr[11], tr[15]));
}
COREDLL void pbrtRotate(float angle, float dx, float dy, float dz) {
	VERIFY_INITIALIZED("Rotate");
	curTransform = curTransform * Rotate(angle, Vector(dx, dy, dz));
}
COREDLL void pbrtScale(float sx, float sy, float sz) {
	VERIFY_INITIALIZED("Scale");
	curTransform = curTransform * Scale(sx, sy, sz);
}
COREDLL void pbrtLookAt(float ex, float ey, float ez, float lx, float ly,
	float lz, float ux, float uy, float uz) {
	VERIFY_INITIALIZED("LookAt");
	curTransform = curTransform * LookAt(Point(ex, ey, ez), Point(lx, ly, lz),
		Vector(ux, uy, uz));
}
COREDLL void pbrtCoordinateSystem(const string &name) {
	VERIFY_INITIALIZED("CoordinateSystem");
	namedCoordinateSystems[name] = curTransform;
}
COREDLL void pbrtCoordSysTransform(const string &name) {
	VERIFY_INITIALIZED("CoordSysTransform");
	if (namedCoordinateSystems.find(name) !=
	    namedCoordinateSystems.end())
		curTransform = namedCoordinateSystems[name];
}
COREDLL void pbrtPixelFilter(const string &name,
                           const ParamSet &params) {
	VERIFY_OPTIONS("PixelFilter");
	renderOptions->FilterName = name;
	renderOptions->FilterParams = params;
}
COREDLL void pbrtFilm(const string &type, const ParamSet &params) {
	VERIFY_OPTIONS("Film");
	renderOptions->FilmParams = params;
	renderOptions->FilmName = type;
}
COREDLL void pbrtSampler(const string &name, const ParamSet &params) {
	VERIFY_OPTIONS("Sampler");
	renderOptions->SamplerName = name;
	renderOptions->SamplerParams = params;
}
COREDLL void pbrtAccelerator(const string &name, const ParamSet &params) {
	VERIFY_OPTIONS("Accelerator");
	renderOptions->AcceleratorName = name;
	renderOptions->AcceleratorParams = params;
}
COREDLL void pbrtSurfaceIntegrator(const string &name, const ParamSet &params) {
	VERIFY_OPTIONS("SurfaceIntegrator");
	renderOptions->SurfIntegratorName = name;
	renderOptions->SurfIntegratorParams = params;
}
COREDLL void pbrtVolumeIntegrator(const string &name, const ParamSet &params) {
	VERIFY_OPTIONS("VolumeIntegrator");
	renderOptions->VolIntegratorName = name;
	renderOptions->VolIntegratorParams = params;
}
COREDLL void pbrtCamera(const string &name,
                       const ParamSet &params) {
	VERIFY_OPTIONS("Camera");
	renderOptions->CameraName = name;
	renderOptions->CameraParams = params;
	renderOptions->WorldToCamera = curTransform;
	namedCoordinateSystems["camera"] =
		curTransform.GetInverse();
}
COREDLL void pbrtSearchPath(const string &path) {
	VERIFY_OPTIONS("SearchPath");
	UpdatePluginPath(path);
	renderOptions->gotSearchPath = true;
}
COREDLL void pbrtWorldBegin() {
	VERIFY_OPTIONS("WorldBegin");
	currentApiState = STATE_WORLD_BLOCK;
	curTransform = Transform();
	namedCoordinateSystems["world"] = curTransform;
}
COREDLL void pbrtAttributeBegin() {
	VERIFY_WORLD("AttributeBegin");
	pushedGraphicsStates.push_back(graphicsState);
	pushedTransforms.push_back(curTransform);
}
COREDLL void pbrtAttributeEnd() {
	VERIFY_WORLD("AttributeEnd");
	if (!pushedGraphicsStates.size()) {
		Error("Unmatched pbrtAttributeEnd() encountered. "
			"Ignoring it.");
		return;
	}
	graphicsState = pushedGraphicsStates.back();
	curTransform = pushedTransforms.back();
	pushedGraphicsStates.pop_back();
	pushedTransforms.pop_back();
}
COREDLL void pbrtTransformBegin() {
	VERIFY_WORLD("TransformBegin");
	pushedTransforms.push_back(curTransform);
}
COREDLL void pbrtTransformEnd() {
	VERIFY_WORLD("TransformEnd");
	if (!pushedTransforms.size()) {
		Error("Unmatched pbrtTransformEnd() encountered. "
			"Ignoring it.");
		return;
	}
	curTransform = pushedTransforms.back();
	pushedTransforms.pop_back();
}
COREDLL void pbrtTexture(const string &name,
                         const string &type,
						 const string &texname,
						 const ParamSet &params) {
	VERIFY_WORLD("Texture");
	TextureParams tp(params, params,
	                 graphicsState.floatTextures,
		             graphicsState.spectrumTextures);
	if (type == "float")  {
		// Create _float_ texture and store in _floatTextures_
		if (graphicsState.floatTextures.find(name) !=
		    graphicsState.floatTextures.end())
			Warning("Texture \"%s\" being redefined", name.c_str());
		Reference<Texture<float> > ft = MakeFloatTexture(texname,
			curTransform, tp);
		if (ft) graphicsState.floatTextures[name] = ft;
	}
	else if (type == "color")  {
		// Create _color_ texture and store in _spectrumTextures_
		if (graphicsState.spectrumTextures.find(name) != graphicsState.spectrumTextures.end())
			Warning("Texture \"%s\" being redefined", name.c_str());
		Reference<Texture<Spectrum> > st = MakeSpectrumTexture(texname,
			curTransform, tp);
		if (st) graphicsState.spectrumTextures[name] = st;
	}
	else
		Error("Texture type \"%s\" unknown.", type.c_str());
}
COREDLL void pbrtMaterial(const string &name, const ParamSet &params) {
	VERIFY_WORLD("Material");
	graphicsState.material = name;
	graphicsState.materialParams = params;
}
COREDLL void pbrtLightSource(const string &name,
                             const ParamSet &params) {
	VERIFY_WORLD("LightSource");
	Light *lt = MakeLight(name, curTransform, params);
	if (lt == NULL)
		Error("pbrtLightSource: light type "
		      "\"%s\" unknown.", name.c_str());
	else
		renderOptions->lights.push_back(lt);
}
COREDLL void pbrtAreaLightSource(const string &name,
                                 const ParamSet &params) {
	VERIFY_WORLD("AreaLightSource");
	graphicsState.areaLight = name;
	graphicsState.areaLightParams = params;
}
COREDLL void pbrtShape(const string &name,
                       const ParamSet &params) {
	VERIFY_WORLD("Shape");
	Reference<Shape> shape = MakeShape(name,
		curTransform, graphicsState.reverseOrientation,
		params);
	if (!shape) return;
	params.ReportUnused();
	// Initialize area light for shape
	AreaLight *area = NULL;
	if (graphicsState.areaLight != "")
		area = MakeAreaLight(graphicsState.areaLight,
		curTransform, graphicsState.areaLightParams, shape);
	// Initialize material for shape
	TextureParams mp(params,
	                 graphicsState.materialParams,
					 graphicsState.floatTextures,
					 graphicsState.spectrumTextures);
	Reference<Texture<float> > bump = NULL;
	Reference<Material> mtl =
		MakeMaterial(graphicsState.material,
		             curTransform, mp);
	if (!mtl)
		mtl = MakeMaterial("matte", curTransform, mp);
	if (!mtl)
		Severe("Unable to create \"matte\" material?!");
	// Create primitive and add to scene or current instance
	Reference<Primitive> prim =
		new GeometricPrimitive(shape, mtl, area);
	if (renderOptions->currentInstance) {
		if (area)
			Warning("Area lights not supported "
			        "with object instancing");
		renderOptions->currentInstance->push_back(prim);
	}
	else {
		renderOptions->primitives.push_back(prim);
		if (area != NULL) {
			// Add area light for primitive to light vector
			renderOptions->lights.push_back(area);
		}
	}
}
COREDLL void pbrtReverseOrientation() {
	VERIFY_WORLD("ReverseOrientation");
	graphicsState.reverseOrientation =
		!graphicsState.reverseOrientation;
}
COREDLL void pbrtVolume(const string &name,
                        const ParamSet &params) {
	VERIFY_WORLD("Volume");
	VolumeRegion *vr = MakeVolumeRegion(name,
		curTransform, params);
	if (vr) renderOptions->volumeRegions.push_back(vr);
}
COREDLL void pbrtObjectBegin(const string &name) {
	VERIFY_WORLD("ObjectBegin");
	pbrtAttributeBegin();
	if (renderOptions->currentInstance)
		Error("ObjectBegin called inside "
		      "of instance definition");
	renderOptions->instances[name] =
		vector<Reference<Primitive> >();
	renderOptions->currentInstance =
		&renderOptions->instances[name];
}
COREDLL void pbrtObjectEnd() {
	VERIFY_WORLD("ObjectEnd");
	if (!renderOptions->currentInstance)
		Error("ObjectEnd called outside "
		      "of instance definition");
	renderOptions->currentInstance = NULL;
	pbrtAttributeEnd();
}
COREDLL void pbrtObjectInstance(const string &name) {
	VERIFY_WORLD("ObjectInstance");
	// Object instance error checking
	if (renderOptions->currentInstance) {
		Error("ObjectInstance can't be called inside instance definition");
		return;
	}
	if (renderOptions->instances.find(name) == renderOptions->instances.end()) {
		Error("Unable to find instance named \"%s\"", name.c_str());
		return;
	}
	vector<Reference<Primitive> > &in =
		renderOptions->instances[name];
	if (in.size() == 0) return;
	if (in.size() > 1 || !in[0]->CanIntersect()) {
		// Refine instance _Primitive_s and create aggregate
		Reference<Primitive> accel =
			MakeAccelerator(renderOptions->AcceleratorName,
			               in, renderOptions->AcceleratorParams);
		if (!accel)
			accel = MakeAccelerator("kdtree", in, ParamSet());
		if (!accel)
			Severe("Unable to find \"kdtree\" accelerator");
		in.erase(in.begin(), in.end());
		in.push_back(accel);
	}
	Reference<Primitive> prim = new InstancePrimitive(in[0],
		curTransform);
	renderOptions->primitives.push_back(prim);
}
COREDLL void pbrtWorldEnd() {
	VERIFY_WORLD("WorldEnd");
	// Ensure the search path was set
	if (!renderOptions->gotSearchPath)
		Severe("PBRT_SEARCHPATH environment variable "
		                  "wasn't set and a plug-in\n"
			              "search path wasn't given in the "
						  "input (with the SearchPath "
						  "directive).\n");
	// Ensure there are no pushed graphics states
	while (pushedGraphicsStates.size()) {
		Warning("Missing end to pbrtAttributeBegin()");
		pushedGraphicsStates.pop_back();
		pushedTransforms.pop_back();
	}
	// Create scene and render
	Scene *scene = renderOptions->MakeScene();
	if (scene) scene->Render();
	delete scene;
	// Clean up after rendering
	currentApiState = STATE_OPTIONS_BLOCK;
	StatsPrint(stdout);
	curTransform = Transform();
	namedCoordinateSystems.erase(namedCoordinateSystems.begin(),
		namedCoordinateSystems.end());
}
Scene *RenderOptions::MakeScene() const {
	// Create scene objects from API settings
	Filter *filter = MakeFilter(FilterName, FilterParams);
	Film *film = MakeFilm(FilmName, FilmParams, filter);
	Camera *camera = MakeCamera(CameraName, CameraParams,
		WorldToCamera, film);
	Sampler *sampler = MakeSampler(SamplerName, SamplerParams, film);
	SurfaceIntegrator *surfaceIntegrator = MakeSurfaceIntegrator(SurfIntegratorName,
		SurfIntegratorParams);
	VolumeIntegrator *volumeIntegrator = MakeVolumeIntegrator(VolIntegratorName,
		VolIntegratorParams);
	Primitive *accelerator = MakeAccelerator(AcceleratorName,
		primitives, AcceleratorParams);
	if (!accelerator) {
		ParamSet ps;
		accelerator = MakeAccelerator("kdtree", primitives, ps);
	}
	if (!accelerator)
		Severe("Unable to find \"kdtree\" accelerator");
	// Initialize _volumeRegion_ from volume region(s)
	VolumeRegion *volumeRegion;
	if (volumeRegions.size() == 0)
		volumeRegion = NULL;
	else if (volumeRegions.size() == 1)
		volumeRegion = volumeRegions[0];
	else
		volumeRegion = new AggregateVolume(volumeRegions);
	// Make sure all plugins initialized properly
	if (!camera || !sampler || !film || !accelerator ||
		!filter || !surfaceIntegrator || !volumeIntegrator) {
		Severe("Unable to create scene due "
		       "to missing plug-ins");
		return NULL;
	}
	Scene *ret = new Scene(camera,
	    surfaceIntegrator, volumeIntegrator,
		sampler, accelerator, lights, volumeRegion);
	// Erase primitives, lights, and volume regions from _RenderOptions_
	primitives.erase(primitives.begin(),
	                primitives.end());
	lights.erase(lights.begin(),
	            lights.end());
	volumeRegions.erase(volumeRegions.begin(),
	                  volumeRegions.end());
	return ret;
}
