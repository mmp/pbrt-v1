
/*
 * pbrt source code Copyright(c) 1998-2005 Matt Pharr and Greg Humphreys
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 * (See file License.txt for complete license)
 */

#ifndef PBRT_MATERIAL_H
#define PBRT_MATERIAL_H
// material.h*
#include "pbrt.h"
#include "primitive.h"
#include "texture.h"
#include "color.h"
#include "reflection.h"
// Material Class Declarations
class COREDLL Material : public ReferenceCounted {
public:
	// Material Interface
	virtual BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading) const = 0;
	virtual ~Material();
	static void Bump(Reference<Texture<float> > d, const DifferentialGeometry &dgGeom,
		const DifferentialGeometry &dgShading, DifferentialGeometry *dgBump);
};
#endif // PBRT_MATERIAL_H
