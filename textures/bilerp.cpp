
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// bilerp.cpp*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"
// BilerpTexture Declarations
template <class T>
class BilerpTexture : public Texture<T> {
public:
	// BilerpTexture Public Methods
	BilerpTexture(TextureMapping2D *m,
				  const T &t00, const T &t01,
			      const T &t10, const T &t11) {
		mapping = m;
		v00 = t00;
		v01 = t01;
		v10 = t10;
		v11 = t11;
	}
	~BilerpTexture() {
		delete mapping;
	}
	T Evaluate(const DifferentialGeometry &dg) const {
		float s, t, dsdx, dtdx, dsdy, dtdy;
		mapping->Map(dg, &s, &t, &dsdx, &dtdx, &dsdy, &dtdy);
		return (1-s)*(1-t) * v00 +
		       (1-s)*(  t) * v01 +
			   (  s)*(1-t) * v10 +
			   (  s)*(  t) * v11;
	}
private:
	// BilerpTexture Private Data
	TextureMapping2D *mapping;
	T v00, v01, v10, v11;
};
// BilerpTexture Method Definitions
extern "C" DLLEXPORT Texture<float> * CreateFloatTexture(const Transform &tex2world,
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
	return new BilerpTexture<float>(map,
		tp.FindFloat("v00", 0.f), tp.FindFloat("v01", 1.f),
		tp.FindFloat("v10", 0.f), tp.FindFloat("v11", 1.f));
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
	return new BilerpTexture<Spectrum>(map,
		tp.FindSpectrum("v00", 0.f), tp.FindSpectrum("v01", 1.f),
		tp.FindSpectrum("v10", 0.f), tp.FindSpectrum("v11", 1.f));
}
