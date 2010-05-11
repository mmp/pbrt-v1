
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

// volumegrid.cpp*
#include "volume.h"
// VolumeGrid Declarations
class VolumeGrid : public DensityRegion {
public:
	// VolumeGrid Public Methods
	VolumeGrid(const Spectrum &sa, const Spectrum &ss, float gg,
	 		const Spectrum &emit, const BBox &e, const Transform &v2w,
			int nx, int ny, int nz, const float *d);
	~VolumeGrid() { delete[] density; }
	BBox WorldBound() const { return WorldToVolume.GetInverse()(extent); }
	bool IntersectP(const Ray &r, float *t0, float *t1) const {
		Ray ray = WorldToVolume(r);
		return extent.IntersectP(ray, t0, t1);
	}
	float Density(const Point &Pobj) const;
	float D(int x, int y, int z) const {
		x = Clamp(x, 0, nx-1);
		y = Clamp(y, 0, ny-1);
		z = Clamp(z, 0, nz-1);
		return density[z*nx*ny + y*nx + x];
	}
private:
	// VolumeGrid Private Data
	float *density;
	const int nx, ny, nz;
	const BBox extent;
};
// VolumeGrid Method Definitions
VolumeGrid::VolumeGrid(const Spectrum &sa,
		const Spectrum &ss, float gg,
 		const Spectrum &emit, const BBox &e,
		const Transform &v2w,
		int x, int y, int z, const float *d)
	: DensityRegion(sa, ss, gg, emit, v2w),
	nx(x), ny(y), nz(z), extent(e) {
	density = new float[nx*ny*nz];
	memcpy(density, d, nx*ny*nz*sizeof(float));
}
float VolumeGrid::Density(const Point &Pobj) const {
	if (!extent.Inside(Pobj)) return 0;
	// Compute voxel coordinates and offsets for _Pobj_
	float voxx = (Pobj.x - extent.pMin.x) /
		(extent.pMax.x - extent.pMin.x) * nx - .5f;
	float voxy = (Pobj.y - extent.pMin.y) /
		(extent.pMax.y - extent.pMin.y) * ny - .5f;
	float voxz = (Pobj.z - extent.pMin.z) /
		(extent.pMax.z - extent.pMin.z) * nz - .5f;
	int vx = Floor2Int(voxx);
	int vy = Floor2Int(voxy);
	int vz = Floor2Int(voxz);
	float dx = voxx - vx, dy = voxy - vy, dz = voxz - vz;
	// Trilinearly interpolate density values to compute local density
	float d00 = Lerp(dx, D(vx, vy, vz),     D(vx+1, vy, vz));
	float d10 = Lerp(dx, D(vx, vy+1, vz),  D(vx+1, vy+1, vz));
	float d01 = Lerp(dx, D(vx, vy, vz+1),  D(vx+1, vy, vz+1));
	float d11 = Lerp(dx, D(vx, vy+1, vz+1),D(vx+1, vy+1, vz+1));
	float d0 = Lerp(dy, d00, d10);
	float d1 = Lerp(dy, d01, d11);
	return Lerp(dz, d0, d1);
}
extern "C" DLLEXPORT VolumeRegion *CreateVolumeRegion(const Transform &volume2world,
		const ParamSet &params) {
	// Initialize common volume region parameters
	Spectrum sigma_a = params.FindOneSpectrum("sigma_a", 0.);
	Spectrum sigma_s = params.FindOneSpectrum("sigma_s", 0.);
	float g = params.FindOneFloat("g", 0.);
	Spectrum Le = params.FindOneSpectrum("Le", 0.);
	Point p0 = params.FindOnePoint("p0", Point(0,0,0));
	Point p1 = params.FindOnePoint("p1", Point(1,1,1));
	int nitems;
	const float *data = params.FindFloat("density", &nitems);
	if (!data) {
		Error("No \"density\" values provided for volume grid?");
		return NULL;
	}
	int nx = params.FindOneInt("nx", 1);
	int ny = params.FindOneInt("ny", 1);
	int nz = params.FindOneInt("nz", 1);
	if (nitems != nx*ny*nz) {
		Error("VolumeGrid has %d density values but nx*ny*nz = %d",
			nitems, nx*ny*nz);
		return NULL;
	}
	return new VolumeGrid(sigma_a, sigma_s, g, Le, BBox(p0, p1),
		volume2world, nx, ny, nz, data);
}
