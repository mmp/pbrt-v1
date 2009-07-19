
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

#ifndef PBRT_DYNLOAD_H
#define PBRT_DYNLOAD_H
// dynload.h*
#include "pbrt.h"
// Runtime Loading Declarations
COREDLL void UpdatePluginPath(const string &newpath);
COREDLL Reference<Shape> MakeShape(const string &name,
	const Transform &object2world, bool reverseOrientation, const ParamSet &paramSet);
COREDLL Reference<Material> MakeMaterial(const string &name,
	const Transform &mtl2world, const TextureParams &mp);
COREDLL Reference<Texture<float> > MakeFloatTexture(const string &name,
	const Transform &tex2world, const TextureParams &tp);
COREDLL Reference<Texture<Spectrum> > MakeSpectrumTexture(const string &name,
	const Transform &tex2world, const TextureParams &tp);
COREDLL Light *MakeLight(const string &name,
	const Transform &light2world, const ParamSet &paramSet);
COREDLL AreaLight *MakeAreaLight(const string &name,
	const Transform &light2world,
	const ParamSet &paramSet, const Reference<Shape> &shape);
COREDLL VolumeRegion *MakeVolumeRegion(const string &name,
	const Transform &light2world, const ParamSet &paramSet);
COREDLL SurfaceIntegrator *MakeSurfaceIntegrator(const string &name,
		const ParamSet &paramSet);
COREDLL VolumeIntegrator *MakeVolumeIntegrator(const string &name,
		const ParamSet &paramSet);
COREDLL Primitive *MakeAccelerator(const string &name,
		const vector<Reference<Primitive> > &prims,
		const ParamSet &paramSet);
COREDLL Camera *MakeCamera(const string &name,
	const ParamSet &paramSet, const Transform &world2cam, Film *film);
COREDLL Sampler *MakeSampler(const string &name,
	const ParamSet &paramSet, const Film *film);
COREDLL Filter *MakeFilter(const string &name,
	const ParamSet &paramSet);
COREDLL ToneMap *MakeToneMap(const string &name,
	const ParamSet &paramSet);
COREDLL Film *MakeFilm(const string &name,
	const ParamSet &paramSet, Filter *filt);
#endif // PBRT_DYNLOAD_H
