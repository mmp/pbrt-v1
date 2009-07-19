
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// windy.cpp*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"
// WindyTexture Declarations
template <class T> class WindyTexture : public Texture<T> {
public:
	// WindyTexture Public Methods
	~WindyTexture() {
		delete mapping;
	}
	WindyTexture(TextureMapping3D *map) {
		mapping = map;
	}
	T Evaluate(const DifferentialGeometry &dg) const {
		Vector dpdx, dpdy;
		Point P = mapping->Map(dg, &dpdx, &dpdy);
		float windStrength =
			FBm(.1f * P, .1f * dpdx, .1f * dpdy, .5f, 3);
		float waveHeight =
			FBm(P, dpdx, dpdy, .5f, 6);
		return fabsf(windStrength) * waveHeight;
	}
private:
	// WindyTexture Private Data
	TextureMapping3D *mapping;
};
// WindyTexture Method Definitions
extern "C" DLLEXPORT Texture<float> * CreateFloatTexture(const Transform &tex2world,
		const TextureParams &tp) {
	// Initialize 3D texture mapping _map_ from _tp_
	TextureMapping3D *map = new IdentityMapping3D(tex2world);
	return new WindyTexture<float>(map);
}

extern "C" DLLEXPORT Texture<Spectrum> * CreateSpectrumTexture(const Transform &tex2world,
		const TextureParams &tp) {
	// Initialize 3D texture mapping _map_ from _tp_
	TextureMapping3D *map = new IdentityMapping3D(tex2world);
	return new WindyTexture<Spectrum>(map);
}
