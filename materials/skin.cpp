
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

// skin.cpp*
#include "pbrt.h"
#include "material.h"
// Skin Class Declarations
class Skin : public Material {
public:
	Skin(Reference<Texture<float> > bump) : bumpMap(bump) { }
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const;
	Reference<Texture<float> > bumpMap;
};

// Skin Method Definitions
BSDF *Skin::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading) const {
	// Declare skin coefficients
	static float diffuse[3] = {  0.428425f,  0.301341f,  0.331054f};
	static float xy0[3] =     { -1.131747f, -1.016939f, -0.966018f};
	static float z0[3] =      { -1.209182f, -1.462488f, -1.222419f};
	static float e0[3] =      {  6.421658f,  3.699932f,  3.524889f};
	static float xy1[3] =     { -0.546570f, -0.643533f, -0.638934f};
	static float z1[3] =      {  0.380123f,  0.410559f,  0.437367f};
	static float e1[3] =      {  3.685044f,  4.266495f,  4.539742f};
	static float xy2[3] =     { -0.998888f, -1.020153f, -1.027479f};
	static float z2[3] =      {  0.857998f,  0.703913f,  0.573625f};
	static float e2[3] =      { 64.208486f, 63.919687f, 43.809866f};
	static Spectrum xy[3] = { Spectrum(xy0), Spectrum(xy1), Spectrum(xy2) };
	static Spectrum z[3] = { Spectrum(z0), Spectrum(z1), Spectrum(z2) };
	static Spectrum e[3] = { Spectrum(e0), Spectrum(e1), Spectrum(e2) };
	// Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
	DifferentialGeometry dgs;
	if (bumpMap)
		Bump(bumpMap, dgGeom, dgShading, &dgs);
	else
		dgs = dgShading;
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dgs, dgGeom.nn);
	bsdf->Add(BSDF_ALLOC(Lafortune)(Spectrum(diffuse), 3, xy, xy, z, e,
		BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)));
	return bsdf;
}
extern "C" DLLEXPORT Material * CreateMaterial(const Transform &xform,
		const TextureParams &mp) {
	Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
	return new Skin(bumpMap);
}
