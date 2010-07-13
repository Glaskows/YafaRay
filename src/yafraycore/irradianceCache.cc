/****************************************************************************
 *
 * 		irradianceCache.cc: icTree, icRecord and hemisphere types and
 *      operators implementation.
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

#include <yafraycore/irradianceCache.h>

__BEGIN_YAFRAY

// stratifiedHemisphere METHODS
// ***********************************************************************

vector3d_t stratifiedHemisphere::getDirection(int j, int k) {
	if (j<0 || j>M || k<0 || k>N)
		Y_INFO << "ERROR(stratifiedHemisphere.getDirection): j, k out of bound" << std::endl;
	double s1 = rnd();
	double s2 = rnd();
	float tmp = (j+s1)/M;
	float sinTheta = fSqrt(tmp);
	float phi = M_2PI*(k+s2)/N;
	return vector3d_t(sinTheta * fCos(phi),
					  sinTheta * fSin(phi),
					  fSqrt(1 - tmp));
}

// icREC_t METHODS
// ***********************************************************************
//! number of total sections are nSamples = pi*m^2
icRec_t::icRec_t(int m, float kappa):stratHemi(m), kappa(kappa) {
	sampleRadius = std::numeric_limits<float>::max();
}

icRec_t::icRec_t(int m, float kappa, const surfacePoint_t &sp):surfacePoint_t(sp), stratHemi(m), kappa(kappa) {
	sampleRadius = std::numeric_limits<float>::max();
}

vector3d_t icRec_t::getSampleHemisphere(int j, int k) {
	return changeBasis( stratHemi.getDirection(j, k), NU, NV, Nup );
}

void icRec_t::changeSampleRadius(float newr) {
	// we use minimal distance radius (without clamping for now)
	if (newr < sampleRadius) {
		sampleRadius = newr;
	}
}

float icRec_t::getWeight(const icRec_t &record) const {
	float dot = Nup * record.getNup();
	// if record is pointing to the other side, better not to count his contribution
	if (dot<0.f) {
		return 0.f;
	}
	float epNor = fSqrt(1.f - dot) * NORMALIZATION_TERM;
	float epPos = (P - record.P).length() * 2.f / radius;
	float weight = 1.f - kappa * fmax(epPos, epNor);
	return weight;
}

bound_t icRec_t::getBound() const {
	return bound_t(P - radius, P + radius);
}

void icRec_t::setPixelArea(const diffRay_t &ray) {
	spDifferentials_t diff(*this, ray);
	pArea = fSqrt(diff.projectedPixelArea());
	minProjR = 3.f * pArea;
	maxProjR = 20.f * pArea;
}

void icRec_t::setNup(const vector3d_t &wo) {
	Nup = FACE_FORWARD(Ng, N, wo);
}

void icRec_t::limitRbyGradient() {
	radius = std::min(sampleRadius, std::min(irr.R / transGrad[0].length(),
											 std::min(irr.G / transGrad[1].length(),
													  irr.B / transGrad[2].length())) );
}

void icRec_t::clampR() {
	radius = std::min(std::max(radius, minProjR), maxProjR);
}

void icRec_t::clampGradient() {
	for (int i=0; i<3; i++) {
		transGrad[i] =
				transGrad[i] * std::min(1.f, sampleRadius / minProjR);
	}
}

//
// ***********************************************************************

void icTree_t::add(const icRec_t &rec) {
	lock.writeLock();
	const bound_t &bound = rec.getBound();
	recursiveAdd(&root, treeBound, rec, bound,
				 2*rec.getRadius()*M_SQRT3 ); // 2*r*sqrt(3) = (bound.a - bound.g).length
	lock.unlock();
}

bool icTree_t::icLookup_t::operator()(const point3d_t &p, const icRec_t &sample) { // point p isn't used
	float weight = sample.getWeight(record);
	if (weight > 0.f) { // TODO: see if weight > 0 is correct or should be a small number
		// get weighted irradiance sample = E_i(p) * w_i(p)
		// E_i(p) = E_i + (n_i x n) * drotE_i
		color_t rotGradResult, transGradResult;
		vector3d_t NCross = record.getNup() ^ sample.getNup(); // should be normalized already
		rotGradResult.R = NCross * sample.rotGrad[0];
		rotGradResult.G = NCross * sample.rotGrad[1];
		rotGradResult.B = NCross * sample.rotGrad[2];
		vector3d_t posDif = record.P - sample.P;
		transGradResult.R = posDif * sample.transGrad[0];
		transGradResult.G = posDif * sample.transGrad[1];
		transGradResult.B = posDif * sample.transGrad[2];
		radSamples.push_back( (sample.irr + rotGradResult + transGradResult) * weight );
		totalWeight += weight;
	}
	return true; // when could it be false? example?
}

bool icTree_t::getIrradiance(icRec_t &record) {
	icLookup_t lookupProc(record);
	lookup(record.P, lookupProc); // ads weighted radiance values to vector
	// if there is no good irradiance sample return false
	if (lookupProc.radSamples.size() == 0) {
		return false;
	}
	// calculate Sum(E_i(p) * w_i(p))
	for (unsigned int i=0; i<lookupProc.radSamples.size(); i++) {
		record.irr += lookupProc.radSamples[i];
	}
	record.irr = record.irr / lookupProc.totalWeight; // E(p) = Sum(E_i(p) * w_i(p)) / Sum(w_i(p))
	return true;
}

__END_YAFRAY
