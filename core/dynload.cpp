
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

// dynload.cpp*
#include "dynload.h"
#include "paramset.h"
#include "shape.h"
#include "material.h"
#ifndef WIN32
#ifndef __APPLE__
#include <dlfcn.h>
#else
#include <mach-o/dyld.h>
#endif
#endif
#include <map>
using std::map;
// Runtime Loading Forward Declarations
static string SearchPath(const string &searchpath,  // NOBOOK
                         const string &filename); // NOBOOK
template <class D> D *GetPlugin(const string &name,
		map<string, D *> &loadedPlugins,
		const string &searchPath) {
	if (loadedPlugins.find(name) != loadedPlugins.end())
		return loadedPlugins[name];
	string filename = name;
	// Add platform-specific shared library filename suffix
	#ifdef WIN32
	filename += ".dll";
	#else
	filename += ".so";
	#endif
	string path = SearchPath(searchPath, filename);
	D *plugin = NULL;
	if (path != "")
		loadedPlugins[name] = plugin = new D(path.c_str());
	else
		Error("Unable to find Plugin/DLL for \"%s\"",
			name.c_str());
	return plugin;
}
class ShapePlugin;
class FilterPlugin;
class MaterialPlugin;
class TexturePlugin;
class FilmPlugin;
class LightPlugin;
class AreaLightPlugin;
class VolumeRegionPlugin;
class SurfaceIntegratorPlugin;
class VolumeIntegratorPlugin;
class ToneMapPlugin;
class AcceleratorPlugin;
class CameraPlugin;
class SamplerPlugin;
static string SearchPath(const string &searchpath,
                         const string &filename);
// Runtime Loading Static Data
static string PluginSearchPath;
static map<string, ShapePlugin *> shapePlugins;
static map<string, FilterPlugin *> filterPlugins;
static map<string, MaterialPlugin *> materialPlugins;
static map<string, TexturePlugin *> texturePlugins;
static map<string, FilmPlugin *> filmPlugins;
static map<string, LightPlugin *> lightPlugins;
static map<string, AreaLightPlugin *> arealightPlugins;
static map<string, VolumeRegionPlugin *> volumePlugins;
static map<string, SurfaceIntegratorPlugin *> surf_integratorPlugins;
static map<string, VolumeIntegratorPlugin *> vol_integratorPlugins;
static map<string, ToneMapPlugin *> tonemapPlugins;
static map<string, AcceleratorPlugin *> acceleratorPlugins;
static map<string, CameraPlugin *> cameraPlugins;
static map<string, SamplerPlugin *> samplerPlugins;
// Runtime Loading Local Classes
class Plugin {
public:
	// Plugin Public Methods
	Plugin(const string &fname);
	~Plugin();
	void *GetSymbol(const string &symname);
private:
	// Plugin Private Data
	#if defined(WIN32)
	HMODULE hinstLib;
	#elif defined(__APPLE__)
	NSModule hinstLib;
	#else
	void *hinstLib;
	#endif
	string pluginName;
};
class ShapePlugin : public Plugin {
public:
	// ShapePlugin Constructor
	ShapePlugin(const string &name)
		: Plugin(name) {
		CreateShape =
			(CreateShapeFunc)(GetSymbol("CreateShape"));
	}
	typedef Shape * (*CreateShapeFunc)(const Transform &o2w,
		bool reverseOrientation, const ParamSet &params);
	CreateShapeFunc CreateShape;
};
class MaterialPlugin : public Plugin {
	typedef Material * (*CreateMaterialFunc)(const Transform &shader2world, const TextureParams &);
public:
	MaterialPlugin( const string &name ): Plugin( name )
	{
		CreateMaterial =
			(CreateMaterialFunc)(GetSymbol("CreateMaterial"));
	}
	CreateMaterialFunc CreateMaterial;
};
class TexturePlugin : public Plugin {
	typedef Texture<float> *
		(*CreateFloatFunc)(const Transform &shader2world, const TextureParams &);
	typedef Texture<Spectrum> *
		(*CreateSpectrumFunc)(const Transform &shader2world, const TextureParams &);
public:
	TexturePlugin( const string &name ): Plugin( name )
	{
		CreateFloatTex =
			(CreateFloatFunc)(GetSymbol("CreateFloatTexture"));
		CreateSpectrumTex =
			(CreateSpectrumFunc)(GetSymbol("CreateSpectrumTexture"));
	}
	CreateFloatFunc CreateFloatTex;
	CreateSpectrumFunc CreateSpectrumTex;
};
class LightPlugin : public Plugin {
	typedef Light *(*CreateLightFunc)(
		const Transform &light2world, const ParamSet &params);
public:
	LightPlugin( const string &name ): Plugin( name )
	{
		CreateLight =
			(CreateLightFunc)(GetSymbol("CreateLight"));
	}
	CreateLightFunc CreateLight;
};
class VolumeRegionPlugin : public Plugin {
	typedef VolumeRegion *(*CreateVolumeFunc)(const Transform &volume2world,
		const ParamSet &params);
public:
	VolumeRegionPlugin( const string &name ): Plugin( name )
	{
		CreateVolumeRegion =
			(CreateVolumeFunc)(GetSymbol("CreateVolumeRegion"));
	}
	CreateVolumeFunc CreateVolumeRegion;
};
class AreaLightPlugin : public Plugin {
	typedef AreaLight *(*CreateAreaLightFunc)(
		const Transform &light2world, const ParamSet &params,
		const Reference<Shape> &shape);
public:
	AreaLightPlugin( const string &name ): Plugin( name )
	{
		CreateAreaLight =
			(CreateAreaLightFunc)(GetSymbol("CreateAreaLight"));
	}
	CreateAreaLightFunc CreateAreaLight;
};
class SurfaceIntegratorPlugin : public Plugin {
	typedef SurfaceIntegrator *(*CreateSurfaceIntegratorFunc)(const ParamSet &params);
public:
	SurfaceIntegratorPlugin( const string &name ): Plugin( name )
	{
		CreateSurfaceIntegrator =
			(CreateSurfaceIntegratorFunc)(GetSymbol("CreateSurfaceIntegrator"));
	}
	CreateSurfaceIntegratorFunc CreateSurfaceIntegrator;
};
class VolumeIntegratorPlugin : public Plugin {
	typedef VolumeIntegrator *(*CreateVolumeIntegratorFunc)(const ParamSet &params);
public:
	VolumeIntegratorPlugin( const string &name ): Plugin( name )
	{
		CreateVolumeIntegrator =
			(CreateVolumeIntegratorFunc)(GetSymbol("CreateVolumeIntegrator"));
	}
	CreateVolumeIntegratorFunc CreateVolumeIntegrator;
};
class AcceleratorPlugin : public Plugin {
	typedef Primitive *(*CreateAcceleratorFunc)(const vector<Reference<Primitive> > &prims, const ParamSet &params);
public:
	AcceleratorPlugin( const string &name ): Plugin( name )
	{
		CreateAccelerator =
			(CreateAcceleratorFunc)(GetSymbol("CreateAccelerator"));
	}
	CreateAcceleratorFunc CreateAccelerator;
};
class SamplerPlugin : public Plugin {
	typedef Sampler *(*CreateSamplerFunc)(const ParamSet &params, const Film *film);
public:
	SamplerPlugin( const string &name ): Plugin( name )
	{
		CreateSampler =
			(CreateSamplerFunc)(GetSymbol("CreateSampler"));
	}
	CreateSamplerFunc CreateSampler;
};
class CameraPlugin : public Plugin {
	typedef Camera *(*CreateCameraFunc)(const ParamSet &params, const Transform &world2cam, Film *film);
public:
	CameraPlugin( const string &name ): Plugin( name )
	{
		CreateCamera =
			(CreateCameraFunc)(GetSymbol("CreateCamera"));
	}
	CreateCameraFunc CreateCamera;
};
class FilterPlugin : public Plugin {
	typedef Filter *(*CreateFilterFunc)(const ParamSet &params);
public:
	FilterPlugin( const string &name ): Plugin( name ) {
		CreateFilter =
			(CreateFilterFunc)(GetSymbol("CreateFilter"));
	}
	CreateFilterFunc CreateFilter;
};
class ToneMapPlugin : public Plugin {
	typedef ToneMap *(*CreateToneMapFunc)(const ParamSet &params);
public:
	ToneMapPlugin( const string &name ): Plugin( name ) {
		CreateToneMap =
			(CreateToneMapFunc)(GetSymbol("CreateToneMap"));
	}
	CreateToneMapFunc CreateToneMap;
};
class FilmPlugin : public Plugin {
	typedef Film *(*CreateFilmFunc)(const ParamSet &params, Filter *filter);
public:
	FilmPlugin( const string &name ): Plugin( name ) {
		CreateFilm =
			(CreateFilmFunc)(GetSymbol("CreateFilm"));
	}
	CreateFilmFunc CreateFilm;
};
// Runtime Loading Method Definitions
COREDLL void UpdatePluginPath(const string &newpath) {
	string ret;
	for (u_int i = 0; i < newpath.size(); ++i) {
		if (newpath[i] != '&')
			ret += newpath[i];
		else
			ret += PluginSearchPath;
	}
	PluginSearchPath = ret;
}
COREDLL Reference<Shape> MakeShape(const string &name,
		const Transform &object2world,
		bool reverseOrientation,
		const ParamSet &paramSet) {
	ShapePlugin *plugin =
		GetPlugin<ShapePlugin>(name,
							   shapePlugins,
							   PluginSearchPath);
	if (plugin)
		return plugin->CreateShape(object2world,
		                          reverseOrientation,
								  paramSet);
	return NULL;
}
static string SearchPath(const string &searchpath,
		const string &filename) {
	const char *start = searchpath.c_str();
	const char *end = start;
	while (*start) {
		while (*end && *end != PBRT_PATH_SEP[0])
			++end;
		string component(start, end);

		string fn = component + "/" + filename;
		FILE *f = fopen(fn.c_str(), "r");
		if (f) {
			fclose(f);
			return fn;
		}
		if (*end == PBRT_PATH_SEP[0]) ++end;
		start = end;
	}
	return "";
}
COREDLL Reference<Material> MakeMaterial(const string &name,
		const Transform &mtl2world,
		const TextureParams &mp) {
	MaterialPlugin *plugin = GetPlugin<MaterialPlugin>(name, materialPlugins,
		PluginSearchPath);
	if (plugin) {
		Reference<Material> ret =
			plugin->CreateMaterial(mtl2world, mp);
		mp.ReportUnused();
		return ret;
	}
	return NULL;
}
COREDLL Reference<Texture<float> > MakeFloatTexture(const string &name,
		const Transform &tex2world, const TextureParams &tp) {
	TexturePlugin *plugin = GetPlugin<TexturePlugin>(name, texturePlugins,
		PluginSearchPath);
	if (plugin) {
		Reference<Texture<float> > ret =
			plugin->CreateFloatTex(tex2world, tp);
		tp.ReportUnused();
		return ret;
	}
	return NULL;
}
COREDLL Reference<Texture<Spectrum> > MakeSpectrumTexture(const string &name,
		const Transform &tex2world, const TextureParams &tp) {
	TexturePlugin *plugin = GetPlugin<TexturePlugin>(name, texturePlugins,
		PluginSearchPath);
	if (plugin) {
		Reference<Texture<Spectrum> > ret =
			plugin->CreateSpectrumTex(tex2world, tp);
		tp.ReportUnused();
		return ret;
	}
	return NULL;
}
COREDLL Light *MakeLight(const string &name,
		const Transform &light2world, const ParamSet &paramSet) {
	LightPlugin *plugin = GetPlugin<LightPlugin>(name, lightPlugins, PluginSearchPath);
	if (plugin) {
		Light *ret = plugin->CreateLight(light2world, paramSet);
		paramSet.ReportUnused();
		return ret;
	}
	return NULL;
}
COREDLL AreaLight *MakeAreaLight(const string &name,
		const Transform &light2world, const ParamSet &paramSet,
		const Reference<Shape> &shape) {
	AreaLightPlugin *plugin = GetPlugin<AreaLightPlugin>(name, arealightPlugins,
		PluginSearchPath);
	if (plugin) {
		AreaLight *ret = plugin->CreateAreaLight(light2world, paramSet, shape);
		paramSet.ReportUnused();
		return ret;
	}
	return NULL;
}
COREDLL VolumeRegion *MakeVolumeRegion(const string &name,
		const Transform &volume2world, const ParamSet &paramSet) {
	VolumeRegionPlugin *plugin = GetPlugin<VolumeRegionPlugin>(name, volumePlugins,
		PluginSearchPath);
	if (plugin) {
		VolumeRegion *ret = plugin->CreateVolumeRegion(volume2world, paramSet);
		paramSet.ReportUnused();
		return ret;
	}
	return NULL;
}
COREDLL SurfaceIntegrator *MakeSurfaceIntegrator(const string &name,
		const ParamSet &paramSet) {
	SurfaceIntegratorPlugin *plugin = GetPlugin<SurfaceIntegratorPlugin>(name,
		surf_integratorPlugins, PluginSearchPath);
	if (plugin) {
		SurfaceIntegrator *ret = plugin->CreateSurfaceIntegrator(paramSet);
		paramSet.ReportUnused();
		return ret;
	}
	return NULL;
}
COREDLL VolumeIntegrator *MakeVolumeIntegrator(const string &name,
		const ParamSet &paramSet) {
	VolumeIntegratorPlugin *plugin = GetPlugin<VolumeIntegratorPlugin>(name, vol_integratorPlugins,
		PluginSearchPath);
	if (plugin) {
		VolumeIntegrator *ret = plugin->CreateVolumeIntegrator(paramSet);
		paramSet.ReportUnused();
		return ret;
	}
	return NULL;
}
COREDLL Primitive *MakeAccelerator(const string &name,
		const vector<Reference<Primitive> > &prims,
		const ParamSet &paramSet) {
	AcceleratorPlugin *plugin = GetPlugin<AcceleratorPlugin>(name, acceleratorPlugins,
		PluginSearchPath);
	if (plugin) {
		Primitive *ret = plugin->CreateAccelerator(prims, paramSet);
		paramSet.ReportUnused();
		return ret;
	}
	return NULL;
}
COREDLL Camera *MakeCamera(const string &name,
		const ParamSet &paramSet,
		const Transform &world2cam, Film *film) {
	CameraPlugin *plugin = GetPlugin<CameraPlugin>(name, cameraPlugins, PluginSearchPath);
	if (plugin) {
		Camera *ret = plugin->CreateCamera(paramSet, world2cam, film);
		paramSet.ReportUnused();
		return ret;
	}
	return NULL;
}
COREDLL Sampler *MakeSampler(const string &name,
		const ParamSet &paramSet, const Film *film) {
	SamplerPlugin *plugin = GetPlugin<SamplerPlugin>(name, samplerPlugins, PluginSearchPath);
	if (plugin) {
		Sampler *ret = plugin->CreateSampler(paramSet, film);
		paramSet.ReportUnused();
		return ret;
	}
	return NULL;
}
COREDLL Filter *MakeFilter(const string &name,
	const ParamSet &paramSet) {
	FilterPlugin *plugin = GetPlugin<FilterPlugin>(name, filterPlugins, PluginSearchPath);
	if (plugin) {
		Filter *ret = plugin->CreateFilter(paramSet);
		paramSet.ReportUnused();
		return ret;
	}
	return NULL;
}
COREDLL ToneMap *MakeToneMap(const string &name,
	const ParamSet &paramSet) {
	ToneMapPlugin *plugin = GetPlugin<ToneMapPlugin>(name, tonemapPlugins, PluginSearchPath);
	if (plugin) {
		ToneMap *ret = plugin->CreateToneMap(paramSet);
		paramSet.ReportUnused();
		return ret;
	}
	return NULL;
}
COREDLL Film *MakeFilm(const string &name,
	const ParamSet &paramSet, Filter *filter) {
	FilmPlugin *plugin = GetPlugin<FilmPlugin>(name, filmPlugins, PluginSearchPath);
	if (plugin) {
		Film *ret = plugin->CreateFilm(paramSet, filter);
		paramSet.ReportUnused();
		return ret;
	}
	return NULL;
}
// Plugin Method Definitions
Plugin::Plugin(const string &fname) {
	pluginName = fname;
#if defined(WIN32)
	hinstLib = LoadLibrary(fname.c_str());
	if (!hinstLib)
		Error("Can't open plug-in \"%s\"", fname.c_str());
#elif defined(__APPLE__)
	NSObjectFileImage fileImage;
	NSObjectFileImageReturnCode returnCode =
		NSCreateObjectFileImageFromFile(fname.c_str(), &fileImage);
	if(returnCode == NSObjectFileImageSuccess) {
		hinstLib = NSLinkModule(fileImage,fname.c_str(),
			NSLINKMODULE_OPTION_RETURN_ON_ERROR
		     |  NSLINKMODULE_OPTION_PRIVATE);
		NSDestroyObjectFileImage(fileImage);
		if (!hinstLib) {
			Error("Can't open plug-in \"%s\"", fname.c_str());
		}
	}
	else {
		Error("Can't open plug-in \"%s\"", fname.c_str());
	}
#else
	hinstLib = dlopen(fname.c_str(), RTLD_LAZY);
	if (!hinstLib)
		Error("Can't open plug-in \"%s\" (%s)", fname.c_str(),
			dlerror());
#endif
}
Plugin::~Plugin() {
#if defined(WIN32)
	FreeLibrary(hinstLib);
#elif defined(__APPLE__)
	NSUnLinkModule(hinstLib,0);
#else
	dlclose(hinstLib);
#endif
}
void *Plugin::GetSymbol(const string &symname) {
	void *data;
#if defined(WIN32)
	data = GetProcAddress(hinstLib, symname.c_str());
#elif defined(__APPLE__)
	string apple_lossage = string("_") + symname;
	NSSymbol nssym = NSLookupSymbolInModule(hinstLib,apple_lossage.c_str());
	data = NSAddressOfSymbol(nssym);
#else
	data = dlsym(hinstLib, symname.c_str());
#endif
	if (!data)
		Severe("Couldn't get symbol \"%s\" in Plugin %s.", symname.c_str(), pluginName.c_str());
	return data;
}
