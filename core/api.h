
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

#ifndef PBRT_API_H
#define PBRT_API_H
// api.h*
#include "pbrt.h"
// API Function Declarations
extern COREDLL void pbrtIdentity();
extern COREDLL void pbrtTranslate(float dx, float dy, float dz);
extern COREDLL void pbrtRotate(float angle,
                               float ax,
							   float ay,
							   float az);
extern COREDLL void pbrtScale(float sx,
                              float sy,
							  float sz);
extern COREDLL void pbrtLookAt(float ex,
                               float ey,
							   float ez,
							   float lx,
							   float ly,
							   float lz,
							   float ux,
							   float uy,
							   float uz);
extern COREDLL
	void pbrtConcatTransform(float transform[16]);
extern COREDLL
	void pbrtTransform(float transform[16]);
extern COREDLL void pbrtCoordinateSystem(const string &);
extern COREDLL void pbrtCoordSysTransform(const string &);
extern COREDLL void pbrtPixelFilter(const string &name, const ParamSet &params);
extern COREDLL void pbrtFilm(const string &type,
                            const ParamSet &params);
extern COREDLL void pbrtSampler(const string &name,
                               const ParamSet &params);
extern COREDLL void pbrtAccelerator(const string &name,
	                               const ParamSet &params);
extern COREDLL
	void pbrtSurfaceIntegrator(const string &name,
							  const ParamSet &params);
extern COREDLL
	void pbrtVolumeIntegrator(const string &name,
							 const ParamSet &params);
extern COREDLL void pbrtCamera(const string &, const ParamSet &cameraParams);
extern COREDLL void pbrtSearchPath(const string &path);
extern COREDLL void pbrtWorldBegin();
extern COREDLL void pbrtAttributeBegin();
extern COREDLL void pbrtAttributeEnd();
extern COREDLL void pbrtTransformBegin();
extern COREDLL void pbrtTransformEnd();
extern COREDLL void pbrtTexture(const string &name, const string &type,
	const string &texname, const ParamSet &params);
extern COREDLL void pbrtMaterial(const string &name,
                               const ParamSet &params);
extern COREDLL void pbrtLightSource(const string &name, const ParamSet &params);
extern COREDLL void pbrtAreaLightSource(const string &name, const ParamSet &params);
extern COREDLL void pbrtShape(const string &name, const ParamSet &params);
extern COREDLL void pbrtReverseOrientation();
extern COREDLL void pbrtVolume(const string &name, const ParamSet &params);
extern COREDLL void pbrtObjectBegin(const string &name);
extern COREDLL void pbrtObjectEnd();
extern COREDLL void pbrtObjectInstance(const string &name);
extern COREDLL void pbrtWorldEnd();
#endif // PBRT_API_H
