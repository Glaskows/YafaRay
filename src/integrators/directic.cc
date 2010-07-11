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
#include <sstream>

#include <yafraycore/irradianceCache.h>

__BEGIN_YAFRAY

class YAFRAYPLUGIN_EXPORT directIC_t : public mcIntegrator_t
{
public:
	directIC_t(bool transpShad=false, int shadowDepth=4, int rayDepth=6);
	virtual bool preprocess();
	color_t getRadiance(renderState_t &state, ray_t &dir) const;
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

/*! returns the incoming radiance from dir direction
  */
// in dir we return the length of hit
// ray must have from and dir
color_t directIC_t::getRadiance(renderState_t &state, ray_t &ray) const
{
	color_t result;
	ray.tmin = MIN_RAYDIST;
	ray.tmax = -1.0;
	surfacePoint_t hitpoint;
	if (scene->intersect(ray, hitpoint)) {
		BSDF_t matBSDF;
		hitpoint.material->initBSDF(state, hitpoint, matBSDF);
		vector3d_t wo = -ray.dir;
		if (! (matBSDF & BSDF_EMIT) ) {
			if ( matBSDF & (BSDF_DIFFUSE | BSDF_GLOSSY) ) {
				// Totally diffusive! not taking into account any glossyness
				result = estimateAllDirectLight(state, hitpoint, wo);
			}
		}
	} else {
		if (background) {
			result = (*background)(ray, state, false);
		}
		ray.tmax = std::numeric_limits<float>::max();
	}
	return result;
}

colorA_t directIC_t::integrate(renderState_t &state, diffRay_t &ray) const
{
	color_t col(0.0);
	float alpha = 0.0;
	void *o_udat = state.userdata;
	bool oldIncludeLights = state.includeLights;
	icRec_t icRecord(25, 2.5f); // M, Kappa
	// Shoot ray into scene
	float oldRayLength[icRecord.getM()];
	color_t oldRad[icRecord.getM()];
	if(scene->intersect(ray, icRecord)) // If it hits
	{
		// create new memory on stack for material setup
		unsigned char *newUserData[USER_DATA_SIZE];
		state.userdata = (void *)newUserData;
		const material_t *material = icRecord.material;
		BSDF_t bsdfs;
		material->initBSDF(state, icRecord, bsdfs);

		vector3d_t wo = -ray.dir;
		icRecord.setNup(wo);

		if(state.raylevel == 0) state.includeLights = true;

		// obtain material self emitance
		if(bsdfs & BSDF_EMIT) col += material->emit(state, icRecord, wo);

		if(bsdfs & BSDF_DIFFUSE)
		{
			// obtain direct illumination
			col += estimateAllDirectLight(state, icRecord, wo);
			col += estimateCausticPhotons(state, icRecord, wo);
			if(useAmbientOcclusion) col += sampleAmbientOcclusion(state, icRecord, wo);

			// check for an interpolated result
			if (ray.hasDifferentials) {
				icRecord.setPixelArea(ray);
				if (!irrCache->getIrradiance(icRecord)) {
					// we set the projected pixel area on the surface point
					ray_t sRay; // ray from hitpoint to hemisphere sample direction
					sRay.from = icRecord.P;
					//vector3d_t swi; // incoming radiance direction, useless for lambertian diffuse component
					color_t innerRotValues;
					color_t innerTransValuesU;
					color_t innerTransValuesV;
					color_t radiance;
					for (int k=0; k<icRecord.getN(); k++) {
						//Y_INFO << "K = " << k << "********************"<< std::endl;
						innerRotValues.black();
						innerTransValuesU.black();
						innerTransValuesV.black();
						for (int j=0; j<icRecord.getM(); j++) {
							//Y_INFO << "J = " << j << "-------------------"<< std::endl;
							// Calculate each incoming radiance of hemisphere at point icRecord
							sRay.dir = icRecord.getSampleHemisphere(j, k);
							radiance = getRadiance(state, sRay);
							// note: oldRad[j] and oldRayLength[j] means L_j,k-1 and r_j,k-1 respectively
							//       oldRad[j-1] and oldRayLength[j-1] means L_j-1,k and r_j-1,k respectively
							if (k>0) {
								if (j>0) {
								// cos2(theta_j-)sin(theta_j-) * (L_j,k - L_j-1,k) / min(r_j,k , r_j-1,k)
									innerTransValuesU +=
											(icRecord.stratHemi.getCosThetaMinus(j)*icRecord.stratHemi.getCosThetaMinus(j)) /
											(fmin(sRay.tmax, oldRayLength[j-1])) *
											(radiance - oldRad[j-1]);
								//}
								// cos(theta_j)[cos(theta_j-) - cos(theta_j+)] * (L_j,k - L_j,k-1) / [sin(theta_j,k) * min(r_j,k , r_j-1,k)]
								innerTransValuesV +=
										icRecord.stratHemi.getCosTheta(j) * (icRecord.stratHemi.getCosThetaMinus(j) - icRecord.stratHemi.getCosThetaPlus(j)) /
										(icRecord.stratHemi.getSinTheta(j) * fmin(sRay.tmax, oldRayLength[j])) *
										(radiance - oldRad[j]);
							}
							}
							icRecord.irr += radiance;
							icRecord.changeSampleRadius(sRay.tmax);
							// copy new rays and irradiance values over old ones
							oldRad[j] = radiance;
							oldRayLength[j] = sRay.tmax; // ToDo: check that ray length when no intercept is big
							innerRotValues -= icRecord.stratHemi.getTanTheta(j) * radiance;
							//Y_INFO << "Tan(theta): " << icRecord.stratHemi.getTanTheta(j) << " - Radiance: "<< radiance << "\n";
							//Y_INFO << "innerRotValues: " << innerRotValues << std::endl;
						}
						icRecord.rotGrad[0] += icRecord.stratHemi.getVk(k) * innerRotValues.R;
						icRecord.rotGrad[1] += icRecord.stratHemi.getVk(k) * innerRotValues.G;
						icRecord.rotGrad[2] += icRecord.stratHemi.getVk(k) * innerRotValues.B;

						icRecord.transGrad[0] +=
								( (innerTransValuesU.R * M_2PI / (float)icRecord.getN() ) *
								  icRecord.stratHemi.getUk(k) ) +
								( innerTransValuesV.R *
								  icRecord.stratHemi.getVkMinus(k) );
						icRecord.transGrad[1] +=
								( (innerTransValuesU.G * M_2PI / (float)icRecord.getN() ) *
								  icRecord.stratHemi.getUk(k) ) +
								( innerTransValuesV.G *
								  icRecord.stratHemi.getVkMinus(k) );
						icRecord.transGrad[2] +=
								( (innerTransValuesU.B * M_2PI / (float)icRecord.getN() ) *
								  icRecord.stratHemi.getUk(k) ) +
								( innerTransValuesV.B *
								  icRecord.stratHemi.getVkMinus(k) );
						//Y_INFO << "V(k) = " << icRecord.stratHemi.getVk(k) << std::endl;
						//Y_INFO << "Rotgrad(1) = " << icRecord.rotGrad[0] << "\n";
						//Y_INFO << "Rotgrad(2) = " << icRecord.rotGrad[1] << "\n";
						//Y_INFO << "Rotgrad(3) = " << icRecord.rotGrad[2] << "\n";
					}
					float k = M_PI / ((float)icRecord.getM() * (float)icRecord.getN());
					icRecord.irr = icRecord.irr * k;
					//Y_INFO << "k: " << k << " - icRecord.irr: " << icRecord.irr << std::endl;
					//Y_INFO << "NU: "<< icRecord.NU <<  "NV: " << icRecord.NV << "Nup: " << icRecord.getNup() << std::endl;
					for (int i=0; i<3; i++) {
						//Y_INFO << "Rotgrad(" << i << ") = " << icRecord.rotGrad[i] << "\n";
						icRecord.rotGrad[i] = changeBasis(
								icRecord.rotGrad[i] * k,
								icRecord.NU,
								icRecord.NV,
								icRecord.getNup());
						icRecord.transGrad[i] = changeBasis(
								icRecord.transGrad[i],
								icRecord.NU,
								icRecord.NV,
								icRecord.getNup());
						//Y_INFO << "Rotgrad(" << i << ") = " << icRecord.rotGrad[i] << "\n";
						// limit gradient (too big on corners): watch out! circular dependences
						//icRecord.transGrad[i] = icRecord.transGrad[i] * fmin(1.f, icRecord.sampleRadius / icRecord.minProjR);
					}
					// HEURISTICS
					// limit r_i by gradient
					icRecord.radius = std::min(icRecord.sampleRadius,
											   std::min(icRecord.irr.R / icRecord.transGrad[0].length(),
														std::min(icRecord.irr.G / icRecord.transGrad[1].length(),
																 icRecord.irr.B / icRecord.transGrad[2].length()))
											   );
					// clamp r_i
					icRecord.radius = std::min(std::max(icRecord.radius, icRecord.minProjR), icRecord.maxProjR);
					// limit gradient
					for (int i=0; i<3; i++) {
						icRecord.transGrad[i] =
								icRecord.transGrad[i] * std::min(1.f, icRecord.sampleRadius / icRecord.minProjR);
					}
					// END HEURISTICS
					irrCache->add(icRecord);
				}
				col += icRecord.irr * icRecord.material->eval(state, icRecord, wo, icRecord.getNup(), BSDF_DIFFUSE) * M_1_PI;
				//}
			} else
				Y_INFO << "NO DIFFERENTIALS!!!" << std::endl;
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
	//col.clampRGB01();
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

