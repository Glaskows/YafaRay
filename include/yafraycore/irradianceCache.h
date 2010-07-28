/****************************************************************************
 *
 * 		irradianceCache.h: icTree, icRecord and hemisphere definitions.
 *      This is part of the yafaray package
 *      Copyright (C) 2010  George Laskowsky Ziguilinsky
 *
 *      This library is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public
 *      License as published by the Free Software Foundation; either
 *      version 2.1 of the License, or (at your option) any later version.
 *
 *      This library is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *      Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this library; if not, write to the Free Software
 *      Foundation,Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#ifndef Y_IRRADIANCECACHE_H
#define Y_IRRADIANCECACHE_H

#include <core_api/vector3d.h>
#include <core_api/bound.h>
#include <core_api/color.h>
#include <core_api/surface.h>
#include <core_api/material.h>
#include <core_api/scene.h>
#include <yafraycore/octree.h>
#include <utilities/mathOptimizations.h>
#include <utilities/mcqmc.h>

#include <limits>
#include <ctime>

__BEGIN_YAFRAY

inline vector3d_t changeBasis(const vector3d_t &vec, const vector3d_t &nx, const vector3d_t &ny, const vector3d_t &nz) {
	vector3d_t dir(nx * vec.x);
	dir += ny * vec.y;
	dir += nz * vec.z;
	return dir;
};

//! Stratified hemisphere for getting random direction vectors distributed proportionally to the cosine term.
/*!
  * There is no kind of copy constructor or assignment operator, please, don't do that, just use and discard!
  * TODO: checks for boundries in methods! (0<=j<M, 0<=k<N)
  */
struct stratifiedHemisphere {
	stratifiedHemisphere(int nm):
			M(nm), N(M_PI * M), rnd((unsigned)1/*time(0)*/) {
	}
	// METHODS
	vector3d_t getDirection(int j, int k); //!< get random direction sample from section j,k in local coordinate system
	vector3d_t getVk(int k) const; //!< get vector v_k: base-plane vector in the direction (pi/2, phi_k + pi/2)
	vector3d_t getVkMinus(int k) const; //!< get vector v_k: base-plane vector in the direction (pi/2, phi_k + pi/2)
	vector3d_t getUk(int k) const; //!< get vector u_k: base-plane vector in the direction (pi/2, phi_k)
	float getTanTheta(int j) const; //!< get tan(theta_j)
	float getSinTheta(int j) const;
	float getSinThetaMinus(int j) const;
	float getCosTheta(int j) const;
	float getCosThetaMinus(int j) const;
	float getCosThetaPlus(int j) const;

	// VARIABLES
	int M; //!< number of divisions along theta
	int N; //!< number of divisions along phi
	random_t rnd; //!< random number generator. ToDo: change seeds!
};


//! Record of Irradiance Cache (1 bounce and direct lighting)
struct icRec_t : public surfacePoint_t
{
	icRec_t(int m, float kappa);
	icRec_t(int m, float kappa, const surfacePoint_t &sp);
	// METHODS
	vector3d_t		getSampleHemisphere(int j, int k); //!< compute indirect light with direct lighting of first bounce
	float			getWeight(const icRec_t &record) const;
	inline float	getRadius() const { return rClamp; } //!< return the radius of the sample "action" area
	bound_t			getBound() const; //!< get the bounding box of the sample sphere
	int				getM() const { return stratHemi.M; } //!< return the number of division if hemisphere along theta
	int				getN() const { return stratHemi.N; } //!< return the number of division if hemisphere along phi
	void			changeSampleRadius(float newr); //!< change radius if it is a new minimum
	void			setPixelArea(const diffRay_t &dir); //!< calculates the projected pixel area on the surface position of the sample
	void			setNup(const vector3d_t &wo);
	const vector3d_t& getNup() const { return Nup; }
	bool			inFront(const icRec_t &record) const;
	void			clampRbyGradient();
	void			clampRbyScreenSpace();
	void			clampRbyNeighbor();
	void			clampGradient();
	void			setRNeighbor(float r);
	// VARIABLES
	color_t			irr; //!< cached irradiance
	vector3d_t		rotGrad[3];
	vector3d_t		transGrad[3];
	stratifiedHemisphere stratHemi; //!< sampling hemisphere at point location
	float			r; //!< minimum distance of all rays from hemisphere sampling
	float			rMin; //!< min radius based on screen space (1-3 times projected pixel area)
	float			rMax; //!< max radius based on screen space (20 times projected pixel area)
	float			rNeighbor; //!< saves the neighbor clamped radius
	// ToDo: add rotation and translation gradients
	// ToDo: add distance to surfaces, minimum and maximum spacing threshold (for neighbor clamping)
	// ToDo: adaptative sampling
	static const float NORMALIZATION_TERM = 8.113140441; //!< from T&L weight function normalization term 1/sqrt(1-cos10°) for 10°
private:
	float			rClamp; //!< radius clamped by projected pixel area and gradients
	vector3d_t		Nup; //!< normal vector on the side of the hitting ray
	float			pArea; //!< projected pixel area over the sample
	float			kappa; //!< overall changing accuaracy constant
};


struct icTree_t : public octree_t<icRec_t>
{
	icTree_t(const bound_t &bound, int levels):octree_t<icRec_t>(bound, levels) {}
	//! Get irradiance estimation at point p. Return false if there isn't a cached irradiance sample near.
	bool getIrradiance(icRec_t &record);
	//! Add a new cached irradiance sample
	void add(const icRec_t &record);
	//! Perform neighbor clamping on record
	void neighborClamp(icRec_t &record);
	//! Store the entire tree (with IC records data) into an xml file named fileName
	void saveToXml(const std::string &fileName);
private:
	void recursiveFindNear(octNode_t<icRec_t> *node, const bound_t &nodeBound, const icRec_t &record,
						   std::vector<icRec_t *> &nearRecs, float &minR);
	struct icLookup_t {
		icLookup_t(const icRec_t &rec): record(rec),totalWeight(0.f) {}
		bool operator()(const point3d_t &p, const icRec_t &record);
		std::vector<color_t> radSamples;
		const icRec_t &record;
		float totalWeight;
	};
};

// INLINE FUNCTIONS
inline vector3d_t stratifiedHemisphere::getUk(int k) const {
	float phi = M_2PI * ((float)k + 0.5f) / (float)N;
	return vector3d_t(fCos(phi), fSin(phi), 0.f);
}

inline vector3d_t stratifiedHemisphere::getVk(int k) const {
	float phi = M_2PI * ((float)k + 0.5f) / (float)N;
	return vector3d_t(-fSin(phi), fCos(phi), 0.f);
}

inline vector3d_t stratifiedHemisphere::getVkMinus(int k) const {
	float phi = M_2PI * (float)k / (float)N;
	return vector3d_t(-fSin(phi), fCos(phi), 0.f);
}

inline float stratifiedHemisphere::getTanTheta(int j) const {
	return fSqrt( ((float)j+0.5f) / ((float)M - (float)j - 0.5) );
	//return fTan(fAsin( fSqrt( ((float)j+0.5f) / (float) M ) ));
}

inline float stratifiedHemisphere::getSinTheta(int j) const {
	return fSqrt( ((float)j + 0.5f) / (float)M );
}

inline float stratifiedHemisphere::getSinThetaMinus(int j) const {
	return fSqrt( (float)j / (float)M );
}

inline float stratifiedHemisphere::getCosTheta(int j) const {
	return fSqrt( 1.0f - ((float)j + 0.5f) / (float)M );
}

inline float stratifiedHemisphere::getCosThetaMinus(int j) const {
	return fSqrt( 1.0f - (float)j / (float)M );
}

inline float stratifiedHemisphere::getCosThetaPlus(int j) const {
	return fSqrt( 1.0f - ((float)j + 1.0f) / (float)M );
}
__END_YAFRAY

#endif // Y_IRRADIANCECACHE_H
