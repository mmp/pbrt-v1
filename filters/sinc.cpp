
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

// sinc.cpp*
#include "sampling.h"
#include "paramset.h"
// Sinc Filter Declarations
class LanczosSincFilter : public Filter {
public:
	LanczosSincFilter(float xw,
	                  float yw,
					  float t) : Filter(xw, yw) {
		tau = t;
	}
	float Evaluate(float x, float y) const;
	float Sinc1D(float x) const;
private:
	float tau;
};
// Sinc Filter Method Definitions
float LanczosSincFilter::Evaluate(float x, float y) const{
	return Sinc1D(x * invXWidth) * Sinc1D(y * invYWidth);
}
float LanczosSincFilter::Sinc1D(float x) const {
	x = fabsf(x);
	if (x < 1e-5) return 1.f;
	if (x > 1.)   return 0.f;
	x *= M_PI;
	float sinc = sinf(x * tau) / (x * tau);
	float lanczos = sinf(x) / x;
	return sinc * lanczos;
}
extern "C" DLLEXPORT Filter *CreateFilter(const ParamSet &ps) {
	float xw = ps.FindOneFloat("xwidth", 4.);
	float yw = ps.FindOneFloat("ywidth", 4.);
	float tau = ps.FindOneFloat("tau", 3.f);
	return new LanczosSincFilter(xw, yw, tau);
}
