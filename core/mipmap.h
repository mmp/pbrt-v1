
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

#ifndef PBRT_MIPMAP_H
#define PBRT_MIPMAP_H
// mipmap.h*
#include "pbrt.h"
#include "color.h"
// MIPMap Declarations
typedef enum {
	TEXTURE_REPEAT,
	TEXTURE_BLACK,
	TEXTURE_CLAMP
} ImageWrap;
template <class T> class MIPMap {
public:
	// MIPMap Public Methods
	MIPMap(int xres, int yres, const T *data, bool doTri = false,
		   float maxAniso = 8.f, ImageWrap wrapMode = TEXTURE_REPEAT);
	~MIPMap();
	T Lookup(float s, float t, float width = 0.f) const;
	T Lookup(float s, float t, float ds0, float dt0,
		float ds1, float dt1) const;
private:
	// MIPMap Private Methods
	struct ResampleWeight;
	ResampleWeight *resampleWeights(int oldres, int newres) {
		Assert(newres >= oldres);
		ResampleWeight *wt = new ResampleWeight[newres];
		float filterwidth = 2.f;
		for (int i = 0; i < newres; ++i) {
			// Compute image resampling weights for _i_th texel
			float center = (i + .5f) * oldres / newres;
			wt[i].firstTexel = Floor2Int((center - filterwidth) + 0.5f);
			for (int j = 0; j < 4; ++j) {
				float pos = wt[i].firstTexel + j + .5f;
				wt[i].weight[j] = Lanczos((pos - center) / filterwidth);
			}
			// Normalize filter weights for texel resampling
			float invSumWts = 1.f / (wt[i].weight[0] + wt[i].weight[1] +
				wt[i].weight[2] + wt[i].weight[3]);
			for (int j = 0; j < 4; ++j)
				wt[i].weight[j] *= invSumWts;
		}
		return wt;
	}
	float clamp(float v) { return Clamp(v, 0.f, INFINITY); }
	Spectrum clamp(const Spectrum &v) { return v.Clamp(0.f, INFINITY); }
	const T &texel(int level, int s, int t) const;
	T triangle(int level, float s, float t) const;
	T EWA(float s, float t, float ds0, float dt0, float ds1, float dt1, int level) const;
	// MIPMap Private Data
	bool doTrilinear;
	float maxAnisotropy;
	ImageWrap wrapMode;
	struct ResampleWeight {
		int firstTexel;
		float weight[4];
	};
	BlockedArray<T> **pyramid;
	int nLevels;
	#define WEIGHT_LUT_SIZE 128
	static float *weightLut;
};
// MIPMap Method Definitions
template <class T>
MIPMap<T>::MIPMap(int sres, int tres,
                  const T *img, bool doTri,
				  float maxAniso, ImageWrap wm) {
	doTrilinear = doTri;
	maxAnisotropy = maxAniso;
	wrapMode = wm;
	T *resampledImage = NULL;
	if (!IsPowerOf2(sres) || !IsPowerOf2(tres)) {
		// Resample image to power-of-two resolution
		int sPow2 = RoundUpPow2(sres), tPow2 = RoundUpPow2(tres);
		// Resample image in $s$ direction
		ResampleWeight *sWeights = resampleWeights(sres, sPow2);
		resampledImage = new T[sPow2 * tPow2];
		// Apply _sWeights_ to zoom in $s$ direction
		for (int t = 0; t < tres; ++t) {
			for (int s = 0; s < sPow2; ++s) {
				// Compute texel $(s,t)$ in $s$-zoomed image
				resampledImage[t*sPow2+s] = 0.;
				for (int j = 0; j < 4; ++j) {
					int origS = sWeights[s].firstTexel + j;
					if (wrapMode == TEXTURE_REPEAT)
						origS = Mod(origS, sres);
					else if (wrapMode == TEXTURE_CLAMP)
						origS = Clamp(origS, 0, sres-1);
					if (origS >= 0 && origS < sres)
						resampledImage[t*sPow2+s] += sWeights[s].weight[j] *
							img[t*sres + origS];
				}
			}
		}
		delete[] sWeights;
		// Resample image in $t$ direction
		ResampleWeight *tWeights = resampleWeights(tres, tPow2);
		T *workData = new T[tPow2];
		for (int s = 0; s < sPow2; ++s) {
			for (int t = 0; t < tPow2; ++t) {
				workData[t] = 0.;
				for (int j = 0; j < 4; ++j) {
					int offset = tWeights[t].firstTexel + j;
					if (wrapMode == TEXTURE_REPEAT) offset = Mod(offset, tres);
					else if (wrapMode == TEXTURE_CLAMP) offset = Clamp(offset, 0, tres-1);
					if (offset >= 0 && offset < tres)
						workData[t] += tWeights[t].weight[j] *
							resampledImage[offset*sPow2 + s];
				}
			}
			for (int t = 0; t < tPow2; ++t)
				resampledImage[t*sPow2 + s] = clamp(workData[t]);
		}
		delete[] workData;
		delete[] tWeights;
		img = resampledImage;  // XXX need to delete resampledImage when done...
		sres = sPow2;
		tres = tPow2;
	}
	// Initialize levels of MIPMap from image
	nLevels = 1 + Log2Int(float(max(sres, tres)));
	pyramid = new BlockedArray<T> *[nLevels];
	// Initialize most detailed level of MIPMap
	pyramid[0] = new BlockedArray<T>(sres, tres, img);
	for (int i = 1; i < nLevels; ++i) {
		// Initialize $i$th MIPMap level from $i-1$st level
		int sRes = max(1, pyramid[i-1]->uSize()/2);
		int tRes = max(1, pyramid[i-1]->vSize()/2);
		pyramid[i] = new BlockedArray<T>(sRes, tRes);
		// Filter four texels from finer level of pyramid
		for (int t = 0; t < tRes; ++t)
			for (int s = 0; s < sRes; ++s)
				(*pyramid[i])(s, t) = .25f * (
					texel(i-1, 2*s, 2*t) +
					texel(i-1, 2*s+1, 2*t) +
					texel(i-1, 2*s, 2*t+1) +
					texel(i-1, 2*s+1, 2*t+1));
	}
	if (resampledImage) delete[] resampledImage;
	// Initialize EWA filter weights if needed
	if (!weightLut) {
		weightLut = (float *)AllocAligned(WEIGHT_LUT_SIZE *
			sizeof(float));
		for (int i = 0; i < WEIGHT_LUT_SIZE; ++i) {
			float alpha = 2;
			float r2 = float(i) / float(WEIGHT_LUT_SIZE - 1);
			weightLut[i] = expf(-alpha * r2) - expf(-alpha);
		}
	}
}
template <class T>
const T &MIPMap<T>::texel(int level, int s, int t) const {
	const BlockedArray<T> &l = *pyramid[level];
	// Compute texel $(s,t)$ accounting for boundary conditions
	switch (wrapMode) {
		case TEXTURE_REPEAT:
			s = Mod(s, l.uSize());
			t = Mod(t, l.vSize());
			break;
		case TEXTURE_CLAMP:
			s = Clamp(s, 0, l.uSize() - 1);
			t = Clamp(t, 0, l.vSize() - 1);
			break;
		case TEXTURE_BLACK: {
			static const T black = 0.f;
			if (s < 0 || s >= l.uSize() ||
				t < 0 || t >= l.vSize())
				return black;
			break;
		}
	}
	return l(s, t);
}
template <class T>
MIPMap<T>::~MIPMap() {
	for (int i = 0; i < nLevels; ++i)
		delete pyramid[i];
	delete[] pyramid;
}
template <class T>
T MIPMap<T>::Lookup(float s, float t, float width) const {
	static StatsCounter mipTrilerps("Texture", // NOBOOK
		"Trilinear MIPMap lookups"); // NOBOOK
	++mipTrilerps; // NOBOOK
	// Compute MIPMap level for trilinear filtering
	float level = nLevels - 1 + Log2(max(width, 1e-8f));
	// Perform trilinear interpolation at appropriate MIPMap level
	if (level < 0)
		return triangle(0, s, t);
	else if (level >= nLevels - 1)
		return texel(nLevels-1, 0, 0);
	else {
		int iLevel = Floor2Int(level);
		float delta = level - iLevel;
		return (1.f-delta) * triangle(iLevel, s, t) +
			delta * triangle(iLevel+1, s, t);
	}
}
template <class T>
T MIPMap<T>::triangle(int level, float s, float t) const {
	level = Clamp(level, 0, nLevels-1);
	s = s * pyramid[level]->uSize() - 0.5f;
	t = t * pyramid[level]->vSize() - 0.5f;
	int s0 = Floor2Int(s), t0 = Floor2Int(t);
	float ds = s - s0, dt = t - t0;
	return (1.f-ds)*(1.f-dt) * texel(level, s0, t0) +
		(1.f-ds)*dt * texel(level, s0, t0+1) +
		ds*(1.f-dt) * texel(level, s0+1, t0) +
		ds*dt * texel(level, s0+1, t0+1);
}
template <class T>
T MIPMap<T>::Lookup(float s, float t, float ds0, float dt0,
		float ds1, float dt1) const {
	static StatsCounter ewaLookups("Texture", "EWA filter lookups"); // NOBOOK
	++ewaLookups; // NOBOOK
	if (doTrilinear)
		return Lookup(s, t,
		              2.f * max(max(fabsf(ds0),
					                fabsf(dt0)),
		                        max(fabsf(ds1),
								    fabsf(dt1))));
	// Compute ellipse minor and major axes
	if (ds0*ds0 + dt0*dt0 < ds1*ds1 + dt1*dt1) {
		swap(ds0, ds1);
		swap(dt0, dt1);
	}
	float majorLength = sqrtf(ds0*ds0 + dt0*dt0);
	float minorLength = sqrtf(ds1*ds1 + dt1*dt1);
	// Clamp ellipse eccentricity if too large
	if (minorLength * maxAnisotropy < majorLength && minorLength > 0.f) {
		float scale = majorLength /
		              (minorLength * maxAnisotropy);
		ds1 *= scale;
		dt1 *= scale;
		minorLength *= scale;
	}
	// Choose level of detail for EWA lookup and perform EWA filtering
	float lod = max(0.f, nLevels - 1.f + Log2(minorLength));
	int ilod = Floor2Int(lod);
	float d = lod - ilod;
	return (1.f - d) * EWA(s, t, ds0, dt0, ds1, dt1, ilod) +
		d * EWA(s, t, ds0, dt0, ds1, dt1, ilod+1);
}
template <class T>
T MIPMap<T>::EWA(float s, float t, float ds0, float dt0,
		float ds1, float dt1, int level) const {
	if (level >= nLevels) return texel(nLevels-1, 0, 0);
	// Convert EWA coordinates to appropriate scale for level
	s = s * pyramid[level]->uSize() - 0.5f;
	t = t * pyramid[level]->vSize() - 0.5f;
	ds0 *= pyramid[level]->uSize();
	dt0 *= pyramid[level]->vSize();
	ds1 *= pyramid[level]->uSize();
	dt1 *= pyramid[level]->vSize();
	// Compute ellipse coefficients to bound EWA filter region
	float A = dt0*dt0 + dt1*dt1 + 1;
	float B = -2.f * (ds0*dt0 + ds1*dt1);
	float C = ds0*ds0 + ds1*ds1 + 1;
	float invF = 1.f / (A*C - B*B*0.25f);
	A *= invF;
	B *= invF;
	C *= invF;
	// Compute the ellipse's $(s,t)$ bounding box in texture space
	float det = -B*B + 4.f*A*C;
	float invDet = 1.f / det;
	float uSqrt = sqrtf(det * C), vSqrt = sqrtf(A * det);
	int s0 = Ceil2Int (s - 2.f * invDet * uSqrt);
	int s1 = Floor2Int(s + 2.f * invDet * uSqrt);
	int t0 = Ceil2Int (t - 2.f * invDet * vSqrt);
	int t1 = Floor2Int(t + 2.f * invDet * vSqrt);
	static StatsRatio ewaTexels("Texture", "Texels per EWA lookup"); // NOBOOK
	ewaTexels.Add((1+s1-s0) * (1+t1-t0), 1); // NOBOOK
	// Scan over ellipse bound and compute quadratic equation
	T num(0.);
	float den = 0;
	for (int it = t0; it <= t1; ++it) {
		float tt = it - t;
		for (int is = s0; is <= s1; ++is) {
			float ss = is - s;
			// Compute squared radius and filter texel if inside ellipse
			float r2 = A*ss*ss + B*ss*tt + C*tt*tt;
			if (r2 < 1.) {
				float weight =
					weightLut[min(Float2Int(r2 * WEIGHT_LUT_SIZE),
					WEIGHT_LUT_SIZE-1)];
				num += texel(level, is, it) * weight;
				den += weight;
			}
		}
	}
	return num / den;
}
template <class T> float *MIPMap<T>::weightLut = NULL;
#endif // PBRT_MIPMAP_H
