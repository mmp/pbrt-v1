
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

// checkerboard.cpp*
#include "pbrt.h"
#include "texture.h"
#include "paramset.h"
#include "sampling.h"
#include "shape.h"
// Checkerboard2D Private Types
enum AAMethod { NONE, SUPERSAMPLE, CLOSEDFORM };
// CheckerboardTexture Declarations
template <class T> class Checkerboard2D : public Texture<T> {
public:
	// Checkerboard2D Public Methods
	Checkerboard2D(TextureMapping2D *m,
	               Reference<Texture<T> > c1,
			       Reference<Texture<T> > c2,
				   const string &aa) {
		mapping = m;
		tex1 = c1;
		tex2 = c2;
		// Select anti-aliasing method for _Checkerboard2D_
		if (aa == "none") aaMethod = NONE;
		else if (aa == "supersample") aaMethod = SUPERSAMPLE;
		else if (aa == "closedform") aaMethod = CLOSEDFORM;
		else {
			Warning("Anti-aliasing mode \"%s\" not understood "
				"by Checkerboard2D, defaulting"
				"to \"supersample\"", aa.c_str());
			aaMethod = SUPERSAMPLE;
		}
	}
	~Checkerboard2D() {
		delete mapping;
	}
	T Evaluate(const DifferentialGeometry &dg) const {
		float s, t, dsdx, dtdx, dsdy, dtdy;
		mapping->Map(dg, &s, &t, &dsdx, &dtdx, &dsdy, &dtdy);
		if (aaMethod == CLOSEDFORM) {
			// Compute closed form box-filtered _Checkerboard2D_ value
			// Evaluate single check if filter is entirely inside one of them
			float ds = max(fabsf(dsdx), fabsf(dsdy));
			float dt = max(fabsf(dtdx), fabsf(dtdy));
			float s0 = s - ds, s1 = s + ds;
			float t0 = t - dt, t1 = t + dt;
			if (Floor2Int(s0) == Floor2Int(s1) &&
				Floor2Int(t0) == Floor2Int(t1)) {
				// Point sample _Checkerboard2D_
				if ((Floor2Int(s) + Floor2Int(t)) % 2 == 0)
					return tex1->Evaluate(dg);
				return tex2->Evaluate(dg);
			}
			// Apply box-filter to checkerboard region
			#define BUMPINT(x) \
				(Floor2Int((x)/2) + \
				 2.f * max((x/2)-Floor2Int(x/2) - .5f, 0.f))
			float sint = (BUMPINT(s1) - BUMPINT(s0)) / (2.f * ds);
			float tint = (BUMPINT(t1) - BUMPINT(t0)) / (2.f * dt);
			float area2 = sint + tint - 2.f * sint * tint;
			if (ds > 1.f || dt > 1.f)
				area2 = .5f;
			return (1.f - area2) * tex1->Evaluate(dg) +
				area2 * tex2->Evaluate(dg);
		}
		else if (aaMethod == SUPERSAMPLE) {
			// Supersample _Checkerboard2D_
			#define SQRT_SAMPLES 4
			#define N_SAMPLES (SQRT_SAMPLES * SQRT_SAMPLES)
			float samples[2*N_SAMPLES];
			StratifiedSample2D(samples, SQRT_SAMPLES, SQRT_SAMPLES);
			T value = 0.;
			float filterSum = 0.;
			for (int i = 0; i < N_SAMPLES; ++i) {
				// Compute new differential geometry for supersample location
				float dx = samples[2*i]   - 0.5f;
				float dy = samples[2*i+1] - 0.5f;
				DifferentialGeometry dgs = dg;
				dgs.p += dx * dgs.dpdx + dy * dgs.dpdy;
				dgs.u += dx * dgs.dudx + dy * dgs.dudy;
				dgs.v += dx * dgs.dvdx + dy * dgs.dvdy;
				dgs.dudx /= N_SAMPLES;
				dgs.dudy /= N_SAMPLES;
				dgs.dvdx /= N_SAMPLES;
				dgs.dvdy /= N_SAMPLES;
				// Compute $(s,t)$ for supersample and evaluate sub-texture
				float ss, ts, dsdxs, dtdxs, dsdys, dtdys;
				mapping->Map(dgs, &ss, &ts, &dsdxs, &dtdxs, &dsdys, &dtdys);
				float wt = expf(-2.f * (dx*dx + dy*dy));
				filterSum += wt;
				if ((Floor2Int(ss) + Floor2Int(ts)) % 2 == 0)
					value += wt * tex1->Evaluate(dgs);
				else
					value += wt * tex2->Evaluate(dgs);
			}
			return value / filterSum;
			#undef N_SAMPLES // NOBOOK
		}
		// Point sample _Checkerboard2D_
		if ((Floor2Int(s) + Floor2Int(t)) % 2 == 0)
			return tex1->Evaluate(dg);
		return tex2->Evaluate(dg);
	}
private:
	// Checkerboard2D Private Data
	Reference<Texture<T> > tex1, tex2;
	TextureMapping2D *mapping;
	AAMethod aaMethod;
};
template <class T> class Checkerboard3D : public Texture<T> {
public:
	// Checkerboard3D Public Methods
	Checkerboard3D(TextureMapping3D *m,
	               Reference<Texture<T> > c1,
			       Reference<Texture<T> > c2) {
		mapping = m;
		tex1 = c1;
		tex2 = c2;
	}
	T Evaluate(const DifferentialGeometry &dg) const {
		// Supersample _Checkerboard3D_
		#define N_SAMPLES 4
		float samples[2*N_SAMPLES*N_SAMPLES];
		StratifiedSample2D(samples, N_SAMPLES, N_SAMPLES);
		T value = 0.;
		float filterSum = 0.;
		for (int i = 0; i < N_SAMPLES*N_SAMPLES; ++i) {
			// Compute new differential geometry for supersample location
			float dx = samples[2*i]   - 0.5f;
			float dy = samples[2*i+1] - 0.5f;
			DifferentialGeometry dgs = dg;
			dgs.p += dx * dgs.dpdx + dy * dgs.dpdy;
			dgs.u += dx * dgs.dudx + dy * dgs.dudy;
			dgs.v += dx * dgs.dvdx + dy * dgs.dvdy;
			dgs.dudx /= N_SAMPLES;
			dgs.dudy /= N_SAMPLES;
			dgs.dvdx /= N_SAMPLES;
			dgs.dvdy /= N_SAMPLES;
			// Compute 3D supersample position and evaluate sub-texture
			Vector dPPdx, dPPdy;
			Point PP = mapping->Map(dgs, &dPPdx, &dPPdy);
			float wt = expf(-2.f * (dx*dx + dy*dy));
			filterSum += wt;
			if ((Floor2Int(PP.x) + Floor2Int(PP.y) + Floor2Int(PP.z)) % 2 == 0)
				value += wt * tex1->Evaluate(dgs);
			else
				value += wt * tex2->Evaluate(dgs);
		}
		return value / filterSum;
	}
private:
	// Checkerboard3D Private Data
	Reference<Texture<T> > tex1, tex2;
	TextureMapping3D *mapping;
};
// CheckerboardTexture Method Definitions
extern "C" DLLEXPORT Texture<float> * CreateFloatTexture(const Transform &tex2world,
		const TextureParams &tp) {
	int dim = tp.FindInt("dimension", 2);
	if (dim != 2 && dim != 3) {
		Error("%d dimensional checkerboard texture not supported", dim);
		return NULL;
	}
	Reference<Texture<float> > tex1 = tp.GetFloatTexture("tex1", 1.f);
	Reference<Texture<float> > tex2 = tp.GetFloatTexture("tex2", 0.f);
	if (dim == 2) {
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
		string aamode = tp.FindString("aamode");
		if (aamode == "") aamode = "closedform";
		return new Checkerboard2D<float>(map, tex1, tex2, aamode);
	}
	else {
		// Initialize 3D texture mapping _map_ from _tp_
		TextureMapping3D *map = new IdentityMapping3D(tex2world);
		return new Checkerboard3D<float>(map, tex1, tex2);
	}
}

extern "C" DLLEXPORT Texture<Spectrum> * CreateSpectrumTexture(const Transform &tex2world,
		const TextureParams &tp) {
	int dim = tp.FindInt("dimension", 2);
	if (dim != 2 && dim != 3) {
		Error("%d dimensional checkerboard texture not supported", dim);
		return NULL;
	}
	Reference<Texture<Spectrum> > tex1 = tp.GetSpectrumTexture("tex1", 1.f);
	Reference<Texture<Spectrum> > tex2 = tp.GetSpectrumTexture("tex2", 0.f);
	if (dim == 2) {
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
		string aamode = tp.FindString("aamode");
		if (aamode == "") aamode = "closedform";
		return new Checkerboard2D<Spectrum>(map, tex1, tex2, aamode);
	}
	else {
		// Initialize 3D texture mapping _map_ from _tp_
		TextureMapping3D *map = new IdentityMapping3D(tex2world);
		return new Checkerboard3D<Spectrum>(map, tex1, tex2);
	}
}
