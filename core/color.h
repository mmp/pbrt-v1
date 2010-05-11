
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

#ifndef PBRT_COLOR_H
#define PBRT_COLOR_H
// color.h*
#include "pbrt.h"
// Spectrum Declarations
class COREDLL Spectrum {
public:
	// Spectrum Public Methods
	Spectrum(float v = 0.f) {
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			c[i] = v;
		Assert(!IsNaN());
	}
	Spectrum(float cs[COLOR_SAMPLES]) {
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			c[i] = cs[i];
		Assert(!IsNaN());
	}
        Spectrum(const Spectrum &s2) {
            Assert(!s2.IsNaN());
            for (int i = 0; i < COLOR_SAMPLES; ++i)
                c[i] = s2.c[i];
        }
        Spectrum &operator=(const Spectrum &s2) {
            Assert(!s2.IsNaN());
            for (int i = 0; i < COLOR_SAMPLES; ++i)
                c[i] = s2.c[i];
            return *this;
        }

	friend ostream &operator<<(ostream &, const Spectrum &);
	Spectrum &operator+=(const Spectrum &s2) {
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			c[i] += s2.c[i];
		return *this;
	}
	Spectrum operator+(const Spectrum &s2) const {
		Spectrum ret = *this;
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			ret.c[i] += s2.c[i];
		return ret;
	}
	Spectrum operator-(const Spectrum &s2) const {
		Spectrum ret = *this;
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			ret.c[i] -= s2.c[i];
		return ret;
	}
	Spectrum operator/(const Spectrum &s2) const {
		Spectrum ret = *this;
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			ret.c[i] /= s2.c[i];
		return ret;
	}
	Spectrum operator*(const Spectrum &sp) const {
		Spectrum ret = *this;
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			ret.c[i] *= sp.c[i];
		return ret;
	}
	Spectrum &operator*=(const Spectrum &sp) {
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			c[i] *= sp.c[i];
		return *this;
	}
	Spectrum operator*(float a) const {
		Assert(!isnan(a));
		Spectrum ret = *this;
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			ret.c[i] *= a;
		return ret;
	}
	Spectrum &operator*=(float a) {
		Assert(!isnan(a));
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			c[i] *= a;
		return *this;
	}
	friend inline
	Spectrum operator*(float a, const Spectrum &s) {
		return s * a;
	}
	Spectrum operator/(float a) const {
	        Assert(a != 0.);
		return *this * (1.f / a);
	}
	Spectrum &operator/=(float a) {
		float inv = 1.f / a;
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			c[i] *= inv;
		return *this;
	}
	void AddWeighted(float w, const Spectrum &s) {
		Assert(!isnan(w));
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			c[i] += w * s.c[i];
	}
	bool operator==(const Spectrum &sp) const {
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			if (c[i] != sp.c[i]) return false;
		return true;
	}
	bool operator!=(const Spectrum &sp) const {
		return !(*this == sp);
	}
	bool Black() const {
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			if (c[i] != 0.) return false;
		return true;
	}
	Spectrum Sqrt() const {
		Spectrum ret;
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			ret.c[i] = sqrtf(c[i]);
		Assert(!ret.IsNaN());
		return ret;
	}
	Spectrum Pow(const Spectrum &e) const {
		Spectrum ret;
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			ret.c[i] = c[i] > 0 ? powf(c[i], e.c[i]) : 0.f;
		Assert(!ret.IsNaN());
		return ret;
	}
	Spectrum operator-() const {
		Spectrum ret;
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			ret.c[i] = -c[i];
		return ret;
	}
	friend Spectrum Exp(const Spectrum &s) {
		Spectrum ret;
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			ret.c[i] = expf(s.c[i]);
		Assert(!ret.IsNaN());
		return ret;
	}
	Spectrum Clamp(float low = 0.f,
	               float high = INFINITY) const {
		Spectrum ret;
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			ret.c[i] = ::Clamp(c[i], low, high);
		return ret;
	}
	bool IsNaN() const {
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			if (isnan(c[i])) return true;
		return false;
	}
	void Print(FILE *f) const {
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			fprintf(f, "%f ", c[i]);
	}
	void XYZ(float xyz[3]) const {
		xyz[0] = xyz[1] = xyz[2] = 0.;
		for (int i = 0; i < COLOR_SAMPLES; ++i) {
			xyz[0] += XWeight[i] * c[i];
			xyz[1] += YWeight[i] * c[i];
			xyz[2] += ZWeight[i] * c[i];
		}
	}
	float y() const {
		float v = 0.;
		for (int i = 0; i < COLOR_SAMPLES; ++i)
			v += YWeight[i] * c[i];
		return v;
	}
	bool operator<(const Spectrum &s2) const {
		return y() < s2.y();
	}
	friend class ParamSet;
	
	// Spectrum Public Data
	static const int CIEstart = 360;
	static const int CIEend = 830;
	static const int nCIE = CIEend-CIEstart+1;
	static const float CIE_X[nCIE];
	static const float CIE_Y[nCIE];
	static const float CIE_Z[nCIE];
private:
	// Spectrum Private Data
	float c[COLOR_SAMPLES];
	static float XWeight[COLOR_SAMPLES];
	static float YWeight[COLOR_SAMPLES];
	static float ZWeight[COLOR_SAMPLES];
	friend Spectrum FromXYZ(float x, float y, float z);
};
#endif // PBRT_COLOR_H
