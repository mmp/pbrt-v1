
/*
 * pbrt source code Copyright(c) 1998-2007 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

#ifndef PBRT_TONEMAP_H
#define PBRT_TONEMAP_H
// tonemap.h*
#include "pbrt.h"
#include "film.h"
// ToneMap Declarations
class ToneMap {
public:
	// ToneMap Interface
	virtual ~ToneMap() { }
	virtual void Map(const float *y, int xRes, int yRes,
		float maxDisplayY, float *scale) const = 0;
};
#endif // PBRT_TONEMAP_H
