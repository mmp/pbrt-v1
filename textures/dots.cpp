
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// dots.cpp*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"
// DotsTexture Declarations
template <class T> class DotsTexture : public Texture<T> {
public:
	// DotsTexture Public Methods
	~DotsTexture() {
		delete mapping;
	}
	DotsTexture(TextureMapping2D *m, Reference<Texture<T> > c1,
			Reference<Texture<T> > c2) {
		mapping = m;
		outsideDot = c1;
		insideDot = c2;
	}
	T Evaluate(const DifferentialGeometry &dg) const {
		// Compute cell indices for dots
		float s, t, dsdx, dtdx, dsdy, dtdy;
		mapping->Map(dg, &s, &t, &dsdx, &dtdx, &dsdy, &dtdy);
		int sCell = Floor2Int(s + .5f), tCell = Floor2Int(t + .5f);
		// Return _insideDot_ result if point is inside dot
		if (Noise(sCell+.5f, tCell+.5f) > 0) {
			float radius = .35f;
			float maxShift = 0.5f - radius;
			float sCenter = sCell + maxShift *
				Noise(sCell + 1.5f, tCell + 2.8f);
			float tCenter = tCell + maxShift *
				Noise(sCell + 4.5f, tCell + 9.8f);
			float ds = s - sCenter, dt = t - tCenter;
			if (ds*ds + dt*dt < radius*radius)
				return insideDot->Evaluate(dg);
		}
		return outsideDot->Evaluate(dg);
	}
private:
	// DotsTexture Private Data
	Reference<Texture<T> > outsideDot, insideDot;
	TextureMapping2D *mapping;
};
// DotsTexture Method Definitions
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
	return new DotsTexture<float>(map,
		tp.GetFloatTexture("inside", 1.f),
		tp.GetFloatTexture("outside", 0.f));
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
	return new DotsTexture<Spectrum>(map,
		tp.GetSpectrumTexture("inside", 1.f),
		tp.GetSpectrumTexture("outside", 0.f));
}
