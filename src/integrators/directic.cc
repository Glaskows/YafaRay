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

//#include <core_api/qtfilm.h>

#include <sstream>

__BEGIN_YAFRAY

//! Stratified hemisphere for getting random direction vectors distributed proportionally to the cosine term.
struct stratifiedHemisphere {
	stratifiedHemisphere(int nm, int nn):M(nm),N(nn) {}
	stratifiedHemisphere(int nm):M(nm) { N = M_PI * M; } //!< we get similar division in both direction with N = pi * M
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
struct icRec_t
{
	icRec_t(int m):stratHemi(m) {}//!< number of total sections are nSamples = pi*m^2
	//color_t		estimateIrradiance(color_t *li); //!< calculate the irradiance by averaging all the incoming radiance
	color_t		getSampleHemisphere(int j, int k, const surfacePoint_t &sp); //!< compute indirect light with direct lighting of first bounce
	float		radius; //!< should be proportional to w
	float		w; //!< weight of the sample
	point3d_t	pos; //!< record location
	normal_t	n; //!< surface normal at pos
	color_t		irr; //!< cached irradiance
	// ToDo: add rotation and translation gradients
	// ToDo: add distance to surfaces, minimum and maximum spacing threshold (for neighbor clamping)
	// ToDo: adaptative sampling
	stratifiedHemisphere stratHemi; //!< sampling hemisphere at point location
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

vector3d_t icRec_t::getSampleHemisphere(int j, int k, const surfacePoint_t &sp) {
	vector3d_t dir; // temporal object to avoid innecesary frequent object creation
	//irr.black();
	//surfacePoint_t hitSp;
	//for (int j=0; j<stratHemi.M; j++) {
	//	for (int k=0; k<stratHemi.N; k++) {
			// get stratified sample
			dir = stratHemi.getSample(j, k);
			// convert to world coordinate system
			dir = sp.U * dir.x + sp.V * dir.y + sp.N * dir.z;
			// obtain sample outgoing radiance in sample direction
	//		diffRay_t lRay(sp.P, dir, MIN_RAYDIST);
	//		if (scene.intersect(lRay, hitSp)) {
				// sum radiance to total irradiance
	//			irr += estimateAllDirectLight(state, hitSp, -lRay.dir);
	//		}
	//	}
	//}
}

class icTree_t : public octree_t<icRec_t>
{
	void add(const icRec_t &rec) //NodeData &dat, const bound_t &bound)
	{
		// bounding box that cover record sphere
		bound_t bound(rec.pos-rec.radius, rec.pos+rec.radius);
		lock.writeLock();
		recursiveAdd(&root, treeBound, rec, bound,
					 (bound.a - bound.g).lengthSqr() );
		lock.unlock();
	}
};

class YAFRAYPLUGIN_EXPORT directIC_t: public mcIntegrator_t
{
public:
	directIC_t(bool transpShad=false, int shadowDepth=4, int rayDepth=6);
	virtual bool preprocess();
	virtual colorA_t integrate(renderState_t &state, diffRay_t &ray) const;
	virtual bool render(imageFilm_t *image);
	static integrator_t* factory(paraMap_t &params, renderEnvironment_t &render);
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
	return success;
}

colorA_t directIC_t::integrate(renderState_t &state, diffRay_t &ray) const
{
	color_t col(0.0);
	float alpha = 0.0;
	surfacePoint_t sp;
	void *o_udat = state.userdata;
	bool oldIncludeLights = state.includeLights;

	// temporal variables for indirect illumination calculus
	stratifiedHemisphere sphere(25);
	vector3d_t irrDir;
	surfacePoint_t irrSp;
	color_t irrC;

	// Shoot ray into scene

	if(scene->intersect(ray, sp)) // If it hits
	{
		unsigned char userdata[USER_DATA_SIZE];
		const material_t *material = sp.material;
		BSDF_t bsdfs;

		state.userdata = (void *) userdata;
		vector3d_t wo = -ray.dir;
		if(state.raylevel == 0) state.includeLights = true;

		material->initBSDF(state, sp, bsdfs);

		if(bsdfs & BSDF_EMIT) col += material->emit(state, sp, wo);

		if(bsdfs & BSDF_DIFFUSE)
		{
			col += estimateAllDirectLight(state, sp, wo);
			col += estimateCausticPhotons(state, sp, wo);
			// if(useAmbientOcclusion) col += sampleAmbientOcclusion(state, sp, wo);
			// col += estimateDiffuseIC(state, sp, wo);
			for (int j=0; j<sphere.M; i++) {
				for (int k=0; k<sphere.N; k++) {
					irrDir = sphere.getSample(j, k, sp);
					diffRay_t lRay(sp.P, dir, MIN_RAYDIST);
					if (scene->intersect(lRay, irrSp)) {
						irrC += estimateAllDirectLight(state, irrSp, -lRay.dir);
					}
				}
			}
			irrC = irrC / sphere.M*N;
			col += irrC;
		}

		recursiveRaytrace(state, ray, bsdfs, sp, wo, col, alpha);

		float m_alpha = material->getAlpha(state, sp, wo);
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

bool directIC_t::render(imageFilm_t *image)
{
    //((mcIntegrator_t *)this)->render(image);
    Y_INFO << "HOLA!!!";
    return true;
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

