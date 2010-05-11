
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

// matte.cpp*
#include "pbrt.h"
#include "material.h"
// Matte Class Declarations
class Matte : public Material {
public:
	// Matte Public Methods
	Matte(Reference<Texture<Spectrum> > kd,
			Reference<Texture<float> > sig,
			Reference<Texture<float> > bump) {
		Kd = kd;
		sigma = sig;
		bumpMap = bump;
	}
	BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
	              const DifferentialGeometry &dgShading) const;
private:
	// Matte Private Data
	Reference<Texture<Spectrum> > Kd;
	Reference<Texture<float> > sigma, bumpMap;
};
// Matte Method Definitions
BSDF *Matte::GetBSDF(const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading) const {
	// Allocate _BSDF_, possibly doing bump-mapping with _bumpMap_
	DifferentialGeometry dgs;
	if (bumpMap)
		Bump(bumpMap, dgGeom, dgShading, &dgs);
	else
		dgs = dgShading;
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dgs, dgGeom.nn);
	// Evaluate textures for _Matte_ material and allocate BRDF
	Spectrum r = Kd->Evaluate(dgs).Clamp();
	float sig = Clamp(sigma->Evaluate(dgs), 0.f, 90.f);
	if (sig == 0.)
		bsdf->Add(BSDF_ALLOC(Lambertian)(r));
	else
		bsdf->Add(BSDF_ALLOC(OrenNayar)(r, sig));
	return bsdf;
	return bsdf;
}
extern "C" DLLEXPORT Material * CreateMaterial(const Transform &xform,
		const TextureParams &mp) {
	Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(1.f));
	Reference<Texture<float> > sigma = mp.GetFloatTexture("sigma", 0.f);
	Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
	return new Matte(Kd, sigma, bumpMap);
}
