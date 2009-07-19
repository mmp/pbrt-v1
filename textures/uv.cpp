
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// uv.cpp*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"
// UVTexture Declarations
class UVTexture : public Texture<Spectrum> {
public:
	// UVTexture Public Methods
	UVTexture(TextureMapping2D *m) {
		mapping = m;
	}
	~UVTexture() {
		delete mapping;
	}
	Spectrum Evaluate(const DifferentialGeometry &dg) const {
		float s, t, dsdx, dtdx, dsdy, dtdy;
		mapping->Map(dg, &s, &t, &dsdx, &dtdx, &dsdy, &dtdy);
		float cs[COLOR_SAMPLES];
		memset(cs, 0, COLOR_SAMPLES * sizeof(float));
		cs[0] = s - Floor2Int(s);
		cs[1] = t - Floor2Int(t);
		return Spectrum(cs);
	}
private:
	TextureMapping2D *mapping;
};
// UVTexture Method Definitions
extern "C" DLLEXPORT Texture<float> * CreateFloatTexture(const Transform &tex2world,
		const TextureParams &tp) {
	return NULL;
}

extern "C" DLLEXPORT Texture<Spectrum> * CreateSpectrumTexture(const Transform &tex2world,
		const TextureParams &tp) {
	// Initialize 2D texture mapping _map_ from _tp_
	TextureMapping2D *map = NULL;
	string type = tp.FindString("mapping");
	if (type == "" || type == "uv") {
		float su = tp.FindFloat("uscale", 1.);
		float sv = tp.FindFloat("vscale", 1.);
		float du = tp.FindFloat("udelta", 0.);
		float dv = tp.FindFloat("vdelta", 0.);
		map = new UVMapping2D(su, sv, du, dv);
	}
	else if (type == "spherical") map = new SphericalMapping2D(tex2world.GetInverse());
	else if (type == "cylindrical") map = new CylindricalMapping2D(tex2world.GetInverse());
	else if (type == "planar")
		map = new PlanarMapping2D(tp.FindVector("v1", Vector(1,0,0)),
			tp.FindVector("v2", Vector(0,1,0)),
			tp.FindFloat("udelta", 0.f), tp.FindFloat("vdelta", 0.f));
	else {
		Error("2D texture mapping \"%s\" unknown", type.c_str());
		map = new UVMapping2D;
	}
	return new UVTexture(map);
}
