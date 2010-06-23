/****************************************************************************
 *		directlight.cc: an integrator for direct lighting only
 *		This is part of the yafaray package
 *		Copyright (C) 2006  Mathias Wein (Lynx)
 *		Copyright (C) 2009  Rodrigo Placencia (DarkTide)
 *
 *		This library is free software; you can redistribute it and/or
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
 */

#include <yafray_config.h>
#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/mcintegrator.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <core_api/vector3d.h>
#include <core_api/bound.h>
#include <yafraycore/octree.h>
#include <utilities/mathOptimizations.h>
#include <utilities/mcqmc.h>

#include <limits>

#include <ctime>

//#include <core_api/qtfilm.h>

#include <sstream>

__BEGIN_YAFRAY

//! Stratified hemisphere for getting random direction vectors distributed proportionally to the cosine term.
struct stratifiedHemisphere {
	stratifiedHemisphere(int nm, int nn):M(nm),N(nn) {}
	stratifiedHemisphere(int nm):M(nm), rnd((unsigned)time(0)) { N = M_PI * M; } //!< we get similar division in both direction with N = pi * M
	vector3d_t getSample(int j, int k); //!< get random direction sample from section j,k in local coordinate system
	//int getNumSamples() { return M*N; }
	int M; //!< number of divisions along theta
	int N; //!< number of divisions along phi
	random_t rnd; //!< random number generator. ToDo: change seeds!
};

vector3d_t stratifiedHemisphere::getSample(int j, int k) {
	// ToDo: change seeds!
	double s1 = rnd();
	double s2 = rnd();
	float theta = (j+s1)/M;
	float sinTheta = fSqrt(theta);
	float phi = M_2PI*(k+s2)/N;
	return vector3d_t(sinTheta * fCos(phi),
					  sinTheta * fSin(phi),
					  fSqrt(1 - theta));
}

//! Record of Irradiance Cache (1 bounce and direct lighting)
struct icRec_t : surfacePoint_t
{
	icRec_t(int m, float kappa):stratHemi(m),kappa(kappa) { radius = std::numeric_limits<float>::max(); }//!< number of total sections are nSamples = pi*m^2
	/*color_t		estimateIrradiance(scene_t scene, color_t *li); //!< calculate the irradiance by averaging all the incoming radiance*/
	float			w; //!< weight of the sample
	color_t			irr; //!< cached irradiance
	vector3d_t		getSampleHemisphere(int j, int k, const vector3d_t &wo); //!< compute indirect light with direct lighting of first bounce
	float			getWeight(const surfacePoint_t &p) const;
	float			getRadius() const { return radius; }
	float			getInvRadius() const { return invRadius; }
	bound_t			getBound() const; //!< get the bounding box of the sample sphere
	void			changeRadius(float newr); //!< change radius if it needs too given the new one
	int				getM() const { return stratHemi.M; } //!< return the number of division if hemisphere along theta
	int				getN() const { return stratHemi.N; } //!< return the number of division if hemisphere along phi
	// ToDo: add rotation and translation gradients
	// ToDo: add distance to surfaces, minimum and maximum spacing threshold (for neighbor clamping)
	// ToDo: adaptative sampling
	static const float NORMALIZATION_TERM = 8.113140441; //!< from T&L weight function normalization term 1/sqrt(1-cos10°) for 10°
private:
	stratifiedHemisphere stratHemi; //!< sampling hemisphere at point location
	float			kappa; //!< overall changing accuaracy constant
	float			radius; //!< should be proportional to w
	float			invRadius; //!< inverse of radius
};

/*color_t icRec_t::estimateIrradiance(color_t *li) {
	color_t result;
	int totalSamples = stratHemi.getNumSamples();
	for (int i=0; i<totalSamples; i++) {
		result += li[i];
	}
	result = result / totalSamples;
	return result;
}*/

vector3d_t icRec_t::getSampleHemisphere(int j, int k, const vector3d_t &wo) {
	vector3d_t dir; // temporal object to avoid innecesary frequent object creation
	//irr.black();
	//surfacePoint_t hitSp;
	//for (int j=0; j<stratHemi.M; j++) {
	//	for (int k=0; k<stratHemi.N; k++) {
			// get stratified sample
			dir = stratHemi.getSample(j, k);
			vector3d_t norm = FACE_FORWARD(Ng, N, wo);
			dir = NU*dir.x + NV*dir.y + norm*dir.z;
			// obtain sample outgoing radiance in sample direction
	//		diffRay_t lRay(sp.P, dir, MIN_RAYDIST);
	//		if (scene.intersect(lRay, hitSp)) {
				// sum radiance to total irradiance
	//			irr += estimateAllDirectLight(state, hitSp, -lRay.dir);
	//		}
	//	}
	//}
	return dir;
}

void icRec_t::changeRadius(float newr) {
	// we use minimal distance radius (without clamping for now)
	if (newr < radius) {
		radius = newr;
		invRadius = 1.f / radius;
	}
}

float icRec_t::getWeight(const surfacePoint_t &sp) const {
	float epPos = 2 * (P - sp.P).length() * invRadius;
	float epNor = fSqrt(1 - N * sp.N) * NORMALIZATION_TERM; // TODO: needs FACE_FORWARD restriction?
	return 1 - kappa * fmax(epPos, epNor);
}

bound_t icRec_t::getBound() const {
	return bound_t(P - radius, P + radius);
}

struct icTree_t : public octree_t<icRec_t>
{
	icTree_t(const bound_t &bound):octree_t<icRec_t>(bound) {}
	//! Get irradiance estimation at point p. Return false if there isn't a cached irradiance sample near.
	bool getIrradiance(const surfacePoint_t &p, color_t &irr);
	//! Add a new cached irradiance sample
	void add(const icRec_t &rec)
	{
		lock.writeLock();
		recursiveAdd(&root, treeBound, rec, rec.getBound(),
					 2*rec.getRadius()*M_SQRT3 ); // 2*r*sqrt(3) = (bound.a - bound.g).length
		lock.unlock();
	}
	//bool isNear(surfacePoint_t &p);
	void setBound(const bound_t &bound) { treeBound = bound; }
	struct icLookup_t {
		icLookup_t(const surfacePoint_t &sp): sp(sp) {}
		std::vector<color_t> radSamples;
		const surfacePoint_t sp;
		bool operator()(const point3d_t &p, const icRec_t &);
		float totalWeight;
	};
};

bool icTree_t::icLookup_t::operator()(const point3d_t &p, const icRec_t &sample) {
	float weight = sample.getWeight(sp);
	if (weight > 0.f) { // TODO: see if weight > 0 is correct or should be a small number
		// get weighted irradiance sample = E_i(p) * w_i(p)
		radSamples.push_back( sample.irr * (weight / (sample.getM() * sample.getN()) ) ); // TODO: optimize division?
		totalWeight += weight;
	}
	return true; // when could it be false? example?
}

bool icTree_t::getIrradiance(const surfacePoint_t &sp, color_t &wIrr) {
	icLookup_t lookupProc(sp);
	lookup(sp.P, lookupProc); // ads weighted radiance values to vector
	// if there is no good irradiance sample return false
	if (lookupProc.radSamples.size() == 0)
		return false;
	// calculate Sum(E_i(p) * w_i(p))
	for (unsigned int i=0; i<lookupProc.radSamples.size(); i++) {
		wIrr += lookupProc.radSamples[i];
	}
	wIrr = wIrr / lookupProc.totalWeight; // E(p) = Sum(E_i(p) * w_i(p)) / Sum(w_i(p))
	return true;
}

class YAFRAYPLUGIN_EXPORT directIC_t: public mcIntegrator_t
{
public:
	directIC_t(bool transpShad=false, int shadowDepth=4, int rayDepth=6);
	virtual bool preprocess();
	virtual colorA_t integrate(renderState_t &state, diffRay_t &ray) const;
	static integrator_t* factory(paraMap_t &params, renderEnvironment_t &render);
	icTree_t *irrCache;
};

directIC_t::directIC_t(bool transpShad, int shadowDepth, int rayDepth)
{
	type = SURFACE;
	causRadius = 0.25;
	causDepth = 10;
	nCausPhotons = 100000;
	nCausSearch = 100;
	trShad = transpShad;
	usePhotonCaustics = false;
	sDepth = shadowDepth;
	rDepth = rayDepth;
	intpb = 0;
	integratorName = "DirectIC";
	integratorShortName = "DIC";
}

bool directIC_t::preprocess()
{
	bool success = true;
	std::stringstream set;
	settings = "";

	if(trShad)
	{
		set << "ShadowDepth: [" << sDepth << "]";
	}
	if(!set.str().empty()) set << "+";
	set << "RayDepth: [" << rDepth << "]";

	background = scene->getBackground();
	lights = scene->lights;

	if(usePhotonCaustics)
	{
		success = createCausticMap();
		if(!set.str().empty()) set << "+";
		set << "Caustics:" << nCausPhotons << " photons. ";
	}

	if(useAmbientOcclusion)
	{
		if(!set.str().empty()) set << "+";
		set << "AO";
	}

	settings = set.str();

	//((qtFilm_t *)imageFilm)->drawTriangle(10, 10, 3);

	// setup cache tree
	irrCache = new icTree_t(scene->getSceneBound());

	return success;
}

colorA_t directIC_t::integrate(renderState_t &state, diffRay_t &ray) const
{
	color_t col(0.0);
	float alpha = 0.0;
	void *o_udat = state.userdata;
	bool oldIncludeLights = state.includeLights;

	icRec_t icRecord(5, 5.f); // kappa of 5.0 is temporal, I don't know if this is a good value (1/a = 1/0.20)
	// Shoot ray into scene
	if(scene->intersect(ray, icRecord)) // If it hits
	{
		// create new memory on stack for material setup
		unsigned char *newUserData[USER_DATA_SIZE];
		state.userdata = (void *)newUserData;
		const material_t *material = icRecord.material;
		BSDF_t bsdfs;
		material->initBSDF(state, icRecord, bsdfs);

		vector3d_t wo = -ray.dir;
		if(state.raylevel == 0) state.includeLights = true;

		// obtain material self emitance
		if(bsdfs & BSDF_EMIT) col += material->emit(state, icRecord, wo);

		if(bsdfs & BSDF_DIFFUSE)
		{
			// obtain direct illumination
			col += estimateAllDirectLight(state, icRecord, wo);
			col += estimateCausticPhotons(state, icRecord, wo);
			if(useAmbientOcclusion) col += sampleAmbientOcclusion(state, icRecord, wo);


			ray_t sRay; // ray from hitpoint to hemisphere sample direction
			sRay.from = icRecord.P;
			surfacePoint_t sp; // surface point hit after 1 bounce
			vector3d_t swi; // incoming radiance direction, useless for lambertian diffuse component
			vector3d_t swo;
			for (int j=0; j<icRecord.getM(); j++) {
				for (int k=0; k<icRecord.getN(); k++) {
					// sample ray setup
					sRay.tmin = MIN_RAYDIST;
					sRay.tmax = -1.0;
					sRay.dir = icRecord.getSampleHemisphere(j, k, wo);
					// Calculate each incoming radiance of hemisphere at point icRecord
					if (scene->intersect(sRay, sp)) {
						// we calculate min radius with new value
						icRecord.changeRadius(sRay.tmax);
						BSDF_t matBSDFs = sp.material->getFlags();
						sp.material->initBSDF(state, sp, matBSDFs);
						swo = -sRay.dir;
						if (! (matBSDFs & BSDF_EMIT) ) {
							if (matBSDFs & (BSDF_DIFFUSE | BSDF_GLOSSY)) {
								icRecord.irr += (estimateAllDirectLight(state, sp, swo) *
												icRecord.material->eval(state, icRecord, swo, swi, BSDF_DIFFUSE));
							}
						}
					} else {
						if (background) icRecord.irr += (*background)(sRay, state, false);
					}
				}
			}
			icRecord.irr = icRecord.irr / (icRecord.getM() * icRecord.getN());
			col += icRecord.irr;
		}

		// Reflective?, Refractive?
		recursiveRaytrace(state, ray, bsdfs, icRecord, wo, col, alpha);

		float m_alpha = material->getAlpha(state, icRecord, wo);
		alpha = m_alpha + (1.f - m_alpha) * alpha;
	}
	else // Nothing hit, return background if any
	{
		if(background) col += (*background)(ray, state, false);
	}

	state.userdata = o_udat;
	state.includeLights = oldIncludeLights;
	return colorA_t(col, alpha);
}

integrator_t* directIC_t::factory(paraMap_t &params, renderEnvironment_t &render)
{
	bool transpShad=false;
	bool caustics=false;
	bool do_AO=false;
	int shadowDepth=5;
	int raydepth=5, cDepth=10;
	int search=100, photons=500000;
	int AO_samples = 32;
	double cRad = 0.25;
	double AO_dist = 1.0;
	color_t AO_col(1.f);

	params.getParam("raydepth", raydepth);
	params.getParam("transpShad", transpShad);
	params.getParam("shadowDepth", shadowDepth);
	params.getParam("caustics", caustics);
	params.getParam("photons", photons);
	params.getParam("caustic_mix", search);
	params.getParam("caustic_depth", cDepth);
	params.getParam("caustic_radius", cRad);
	params.getParam("do_AO", do_AO);
	params.getParam("AO_samples", AO_samples);
	params.getParam("AO_distance", AO_dist);
	params.getParam("AO_color", AO_col);

	directIC_t *inte = new directIC_t(transpShad, shadowDepth, raydepth);
	// caustic settings
	inte->usePhotonCaustics = caustics;
	inte->nCausPhotons = photons;
	inte->nCausSearch = search;
	inte->causDepth = cDepth;
	inte->causRadius = cRad;
	// AO settings
	inte->useAmbientOcclusion = do_AO;
	inte->aoSamples = AO_samples;
	inte->aoDist = (PFLOAT)AO_dist;
	inte->aoCol = AO_col;
	return inte;
}

extern "C"
{

	YAFRAYPLUGIN_EXPORT void registerPlugin(renderEnvironment_t &render)
	{
		render.registerFactory("directIC",directIC_t::factory);
	}

}

__END_YAFRAY

