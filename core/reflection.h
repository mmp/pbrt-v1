
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

#ifndef PBRT_REFLECTION_H
#define PBRT_REFLECTION_H
// reflection.h*
#include "pbrt.h"
#include "geometry.h"
#include "shape.h"
// BSDF Inline Functions
inline float CosTheta(const Vector &w) { return w.z; }
inline float SinTheta(const Vector &w) {
	return sqrtf(max(0.f, 1.f - w.z*w.z));
}
inline float SinTheta2(const Vector &w) {
	return max(0.f, 1.f - CosTheta(w)*CosTheta(w));
}
inline float CosPhi(const Vector &w) {
	float sintheta = SinTheta(w);
	if (sintheta == 0.f) return 1.f;
	return Clamp(w.x / sintheta, -1.f, 1.f);
}
inline float SinPhi(const Vector &w) {
	float sintheta = SinTheta(w);
	if (sintheta == 0.f) return 0.f;
	return Clamp(w.y / sintheta, -1.f, 1.f);
}
inline bool SameHemisphere(const Vector &w,
                          const Vector &wp) {
	return w.z * wp.z > 0.f;
}
// BSDF Declarations
enum BxDFType {
	BSDF_REFLECTION   = 1<<0,
	BSDF_TRANSMISSION = 1<<1,
	BSDF_DIFFUSE      = 1<<2,
	BSDF_GLOSSY       = 1<<3,
	BSDF_SPECULAR     = 1<<4,
	BSDF_ALL_TYPES        = BSDF_DIFFUSE |
	                        BSDF_GLOSSY |
					        BSDF_SPECULAR,
	BSDF_ALL_REFLECTION   = BSDF_REFLECTION |
	                        BSDF_ALL_TYPES,
	BSDF_ALL_TRANSMISSION = BSDF_TRANSMISSION |
	                        BSDF_ALL_TYPES,
	BSDF_ALL              = BSDF_ALL_REFLECTION |
	                        BSDF_ALL_TRANSMISSION
};
class COREDLL BSDF {
public:
	// BSDF Public Methods
	Spectrum Sample_f(const Vector &o, Vector *wi, float u1, float u2,
		float u3, float *pdf, BxDFType flags = BSDF_ALL,
		BxDFType *sampledType = NULL) const;
	Spectrum Sample_f(const Vector &wo, Vector *wi,
		BxDFType flags = BSDF_ALL,
		BxDFType *sampledType = NULL) const;
	float Pdf(const Vector &wo,
	          const Vector &wi,
			  BxDFType flags = BSDF_ALL) const;
	BSDF(const DifferentialGeometry &dgs,
	     const Normal &ngeom,
		 float eta = 1.f);
	inline void Add(BxDF *bxdf);
	int NumComponents() const { return nBxDFs; }
	int NumComponents(BxDFType flags) const;
	bool HasShadingGeometry() const {
		return (nn.x != ng.x || nn.y != ng.y || nn.z != ng.z);
	}
	Vector WorldToLocal(const Vector &v) const {
		return Vector(Dot(v, sn), Dot(v, tn), Dot(v, nn));
	}
	Vector LocalToWorld(const Vector &v) const {
		return Vector(sn.x * v.x + tn.x * v.y + nn.x * v.z,
		              sn.y * v.x + tn.y * v.y + nn.y * v.z,
		              sn.z * v.x + tn.z * v.y + nn.z * v.z);
	}
	Spectrum f(const Vector &woW, const Vector &wiW,
		BxDFType flags = BSDF_ALL) const;
	Spectrum rho(BxDFType flags = BSDF_ALL) const;
	Spectrum rho(const Vector &wo,
	             BxDFType flags = BSDF_ALL) const;
	static void *Alloc(u_int sz) { return arena.Alloc(sz); }
	static void FreeAll() { arena.FreeAll(); }
	// BSDF Public Data
	const DifferentialGeometry dgShading;
	const float eta;
private:
	// BSDF Private Methods
	~BSDF() { }
	friend class NoSuchClass;
	// BSDF Private Data
	Normal nn, ng;
	Vector sn, tn;
	int nBxDFs;
	#define MAX_BxDFS 8
	BxDF * bxdfs[MAX_BxDFS];
	static MemoryArena arena;
};
#define BSDF_ALLOC(T)  new (BSDF::Alloc(sizeof(T))) T
// BxDF Declarations
class COREDLL BxDF {
public:
	// BxDF Interface
	virtual ~BxDF() { }
	BxDF(BxDFType t) : type(t) { }
	bool MatchesFlags(BxDFType flags) const {
		return (type & flags) == type;
	}
	virtual Spectrum f(const Vector &wo,
	                   const Vector &wi) const = 0;
	virtual Spectrum Sample_f(const Vector &wo, Vector *wi,
		float u1, float u2, float *pdf) const;
	virtual Spectrum rho(const Vector &wo,
	                     int nSamples = 16,
		                 float *samples = NULL) const;
	virtual Spectrum rho(int nSamples = 16,
	                     float *samples = NULL) const;
	virtual float Pdf(const Vector &wi, const Vector &wo) const;
	// BxDF Public Data
	const BxDFType type;
};
class COREDLL BRDFToBTDF : public BxDF {
public:
	// BRDFToBTDF Public Methods
	BRDFToBTDF(BxDF *b)
		: BxDF(BxDFType(b->type ^
		               (BSDF_REFLECTION | BSDF_TRANSMISSION))) {
		brdf = b;
	}
	static Vector otherHemisphere(const Vector &w) {
		return Vector(w.x, w.y, -w.z);
	}
	Spectrum rho(const Vector &w, int nSamples,
			float *samples) const {
		return brdf->rho(otherHemisphere(w), nSamples, samples);
	}
	Spectrum rho(int nSamples, float *samples) const {
		return brdf->rho(nSamples, samples);
	}
	Spectrum f(const Vector &wo, const Vector &wi) const;
	Spectrum Sample_f(const Vector &wo, Vector *wi,
		float u1, float u2, float *pdf) const;
	float Pdf(const Vector &wo, const Vector &wi) const;
private:
	BxDF *brdf;
};
class COREDLL Fresnel {
public:
	// Fresnel Interface
	virtual ~Fresnel();
	virtual Spectrum Evaluate(float cosi) const = 0;
};
class COREDLL FresnelConductor : public Fresnel {
public:
	// FresnelConductor Public Methods
	Spectrum Evaluate(float cosi) const;
	FresnelConductor(const Spectrum &e, const Spectrum &kk)
		: eta(e), k(kk) {
	}
private:
	// FresnelConductor Private Data
	Spectrum eta, k;
};
class COREDLL FresnelDielectric : public Fresnel {
public:
	// FresnelDielectric Public Methods
	Spectrum Evaluate(float cosi) const;
	FresnelDielectric(float ei, float et) {
		eta_i = ei;
		eta_t = et;
	}
private:
	// FresnelDielectric Private Data
	float eta_i, eta_t;
};
class COREDLL FresnelNoOp : public Fresnel {
public:
	Spectrum Evaluate(float) const { return Spectrum(1.); }
};
class COREDLL SpecularReflection : public BxDF {
public:
	// SpecularReflection Public Methods
	SpecularReflection(const Spectrum &r, Fresnel *f)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)),
		  R(r), fresnel(f) {
	}
	Spectrum f(const Vector &, const Vector &) const {
		return Spectrum(0.);
	}
	Spectrum Sample_f(const Vector &wo, Vector *wi,
	                  float u1, float u2, float *pdf) const;
	float Pdf(const Vector &wo, const Vector &wi) const {
		return 0.;
	}
private:
	// SpecularReflection Private Data
	Spectrum R;
	Fresnel *fresnel;
};
class COREDLL SpecularTransmission : public BxDF {
public:
	// SpecularTransmission Public Methods
	SpecularTransmission(const Spectrum &t, float ei, float et)
		: BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)),
		  fresnel(ei, et) {
		T = t;
		etai = ei;
		etat = et;
	}
	Spectrum f(const Vector &, const Vector &) const {
		return Spectrum(0.);
	}
	Spectrum Sample_f(const Vector &wo, Vector *wi, float u1, float u2, float *pdf) const;
	float Pdf(const Vector &wo, const Vector &wi) const {
		return 0.;
	}
private:
	// SpecularTransmission Private Data
	Spectrum T;
	float etai, etat;
	FresnelDielectric fresnel;
};
class COREDLL Lambertian : public BxDF {
public:
	// Lambertian Public Methods
	Lambertian(const Spectrum &reflectance)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)),
		  R(reflectance), RoverPI(reflectance * INV_PI) {
	}
	Spectrum f(const Vector &wo, const Vector &wi) const;
	Spectrum rho(const Vector &, int, float *) const {
		return R;
	}
	Spectrum rho(int, float *) const { return R; }
private:
	// Lambertian Private Data
	Spectrum R, RoverPI;
};
class COREDLL OrenNayar : public BxDF {
public:
	// OrenNayar Public Methods
	Spectrum f(const Vector &wo, const Vector &wi) const;
	OrenNayar(const Spectrum &reflectance, float sig)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)),
		  R(reflectance) {
		float sigma = Radians(sig);
		float sigma2 = sigma*sigma;
		A = 1.f - (sigma2 / (2.f * (sigma2 + 0.33f)));
		B = 0.45f * sigma2 / (sigma2 + 0.09f);
	}
private:
	// OrenNayar Private Data
	Spectrum R;
	float A, B;
};
class COREDLL MicrofacetDistribution {
public:
	// MicrofacetDistribution Interface
	virtual ~MicrofacetDistribution() { }
	virtual float D(const Vector &wh) const = 0;
	virtual void Sample_f(const Vector &wo, Vector *wi,
		float u1, float u2, float *pdf) const = 0;
	virtual float Pdf(const Vector &wo,
	                  const Vector &wi) const = 0;
};
class COREDLL Microfacet : public BxDF {
public:
	// Microfacet Public Methods
	Microfacet(const Spectrum &reflectance, Fresnel *f,
		MicrofacetDistribution *d);
	Spectrum f(const Vector &wo, const Vector &wi) const;
	float G(const Vector &wo, const Vector &wi,
			const Vector &wh) const {
		float NdotWh = fabsf(CosTheta(wh));
		float NdotWo = fabsf(CosTheta(wo));
		float NdotWi = fabsf(CosTheta(wi));
		float WOdotWh = AbsDot(wo, wh);
		return min(1.f, min((2.f * NdotWh * NdotWo / WOdotWh),
		                (2.f * NdotWh * NdotWi / WOdotWh)));
	}
	Spectrum Sample_f(const Vector &wo, Vector *wi,
		float u1, float u2, float *pdf) const;
	float Pdf(const Vector &wo, const Vector &wi) const;
private:
	// Microfacet Private Data
	Spectrum R;
	MicrofacetDistribution *distribution;
	Fresnel *fresnel;
};
class COREDLL Blinn : public MicrofacetDistribution {
public:
	Blinn(float e) { if (e > 1000.f || isnan(e)) e = 1000.f; exponent = e; }
	// Blinn Public Methods
	float D(const Vector &wh) const {
		float costhetah = fabsf(CosTheta(wh));
		return (exponent+2) *
		       INV_TWOPI *
			   powf(costhetah, exponent);
	}
	virtual void Sample_f(const Vector &wi, Vector *sampled_f, float u1, float u2, float *pdf) const;
	virtual float Pdf(const Vector &wi, const Vector &wo) const;
private:
	float exponent;
};
class COREDLL Anisotropic : public MicrofacetDistribution {
public:
	// Anisotropic Public Methods
	Anisotropic(float x, float y) { ex = x; ey = y;
		if (ex > 1000.f || isnan(ex)) ex = 1000.f;
		if (ey > 1000.f || isnan(ey)) ey = 1000.f;
	}
	float D(const Vector &wh) const {
		float costhetah = fabsf(CosTheta(wh));
		float d = (1.f - costhetah * costhetah);
		if (d == 0.f) return 0.f;
		float e = (ex * wh.x * wh.x + ey * wh.y * wh.y) / d;
		return sqrtf((ex+2.f)*(ey+2.f)) * INV_TWOPI * powf(costhetah, e);
	}
	void Sample_f(const Vector &wo, Vector *wi, float u1, float u2, float *pdf) const;
	float Pdf(const Vector &wo, const Vector &wi) const;
	void sampleFirstQuadrant(float u1, float u2, float *phi, float *costheta) const;
private:
	float ex, ey;
};
class COREDLL Lafortune : public BxDF {
public:
	// Lafortune Public Methods
	Lafortune(const Spectrum &r, u_int nl,
	          const Spectrum *x, const Spectrum *y, const Spectrum *z,
			  const Spectrum *e, BxDFType t);
	Spectrum f(const Vector &wo, const Vector &wi) const;
private:
	// Lafortune Private Data
	Spectrum R;
	u_int nLobes;
	const Spectrum *x, *y, *z, *exponent;
};
class COREDLL FresnelBlend : public BxDF {
public:
	// FresnelBlend Public Methods
	FresnelBlend(const Spectrum &Rd,
	             const Spectrum &Rs,
				 MicrofacetDistribution *dist);
	Spectrum f(const Vector &wo, const Vector &wi) const;
	Spectrum SchlickFresnel(float costheta) const {
		return
			Rs + powf(1 - costheta, 5.f) * (Spectrum(1.) - Rs);
	}
	Spectrum Sample_f(const Vector &wi, Vector *sampled_f, float u1, float u2, float *pdf) const;
	float Pdf(const Vector &wi, const Vector &wo) const;
private:
	// FresnelBlend Private Data
	Spectrum Rd, Rs;
	MicrofacetDistribution *distribution;
};
// BSDF Inline Method Definitions
inline void BSDF::Add(BxDF *b) {
	Assert(nBxDFs < MAX_BxDFS);
	bxdfs[nBxDFs++] = b;
}
inline int BSDF::NumComponents(BxDFType flags) const {
	int num = 0;
	for (int i = 0; i < nBxDFs; ++i)
		if (bxdfs[i]->MatchesFlags(flags)) ++num;
	return num;
}
#endif // PBRT_REFLECTION_H
