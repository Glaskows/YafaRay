/****************************************************************************
 *		mcintegrator.h: A basic abstract integrator for MC sampling
 *		This is part of the yafray package
 *		Copyght (C) 2010  Rodrigo Placencia (DarkTide)
 *
 *		This library is free software; you can redistribute it and/or
 *		modify it under the terms of the GNU Lesser General Public
 *		License as published by the Free Software Foundation; either
 *		version 2.1 of the License, or (at your option) any later version.
 *
 *		This library is distributed in the hope that it will be useful,
 *		but WITHOUT ANY WARRANTY; without even the implied warranty of
 *		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *		Lesser General Public License for more details.
 *
 *		You should have received a copy of the GNU Lesser General Public
 *		License along with this library; if not, write to the Free Software
 *		Foundation,Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include <core_api/mcintegrator.h>
#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <yafraycore/photon.h>
#include <yafraycore/scr_halton.h>
#include <yafraycore/spectrum.h>
#include <utilities/mcqmc.h>

__BEGIN_YAFRAY

#define allBSDFIntersect (BSDF_GLOSSY | BSDF_DIFFUSE | BSDF_DISPERSIVE | BSDF_REFLECT | BSDF_TRANSMIT);
#define loffsDelta 4567 //just some number to have different sequences per light...and it's a prime even...

inline color_t mcIntegrator_t::estimateAllDirectLight(renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo) const
{
	color_t col;
	unsigned int loffs = 0;
	for(std::vector<light_t *>::const_iterator l=lights.begin(); l!=lights.end(); ++l)
	{
		col += doLightEstimation(state, (*l), sp, wo, loffs);
		loffs++;
	}
	
	return col;
}

inline color_t mcIntegrator_t::estimateOneDirectLight(renderState_t &state, const surfacePoint_t &sp, vector3d_t wo, int n) const
{
	int lightNum = lights.size();
	
	if(lightNum == 0) return color_t(0.f); //??? if you get this far the lights must be >= 1 but, what the hell... :)
	
	Halton hal2(2);
	hal2.setStart(n-1);
	
	int lnum = std::min((int)(hal2.getNext() * (float)lightNum), lightNum - 1);
	
	return doLightEstimation(state, lights[lnum], sp, wo, lnum);
	//return col * nLights;
}

inline color_t mcIntegrator_t::doLightEstimation(renderState_t &state, light_t *light, const surfacePoint_t &sp, const vector3d_t &wo, const unsigned int  &loffs) const
{
	color_t col(0.f);
	bool shadowed;
	unsigned int l_offs = loffs * loffsDelta;
	const material_t *material = sp.material;
	ray_t lightRay;
	lightRay.from = sp.P;
	color_t lcol(0.f), scol;
	float lightPdf;

	// handle lights with delta distribution, e.g. point and directional lights
	if( light->diracLight() )
	{
		if( light->illuminate(sp, lcol, lightRay) )
		{
			// ...shadowed...
			lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
			shadowed = (trShad) ? scene->isShadowed(state, lightRay, sDepth, scol) : scene->isShadowed(state, lightRay);
			if(!shadowed)
			{
				if(trShad) lcol *= scol;
				color_t surfCol = material->eval(state, sp, wo, lightRay.dir, BSDF_ALL);
				color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay);
				col += surfCol * lcol * std::fabs(sp.N*lightRay.dir) * transmitCol;
			}
		}
	}
	else // area light and suchlike
	{
		Halton hal2(2);
		Halton hal3(3);
		int n = light->nSamples();
		if(state.rayDivision > 1) n = std::max(1, n/state.rayDivision);
		float invNS = 1.f / (float)n;
		unsigned int offs = n * state.pixelSample + state.samplingOffs + l_offs;
		bool canIntersect=light->canIntersect();
		color_t ccol(0.0);
		lSample_t ls;

		hal2.setStart(offs-1);
		hal3.setStart(offs-1);

		for(int i=0; i<n; ++i)
		{
			// ...get sample val...
			ls.s1 = hal2.getNext();
			ls.s2 = hal3.getNext();
			
			if( light->illumSample (sp, ls, lightRay) )
			{
				// ...shadowed...
				lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is bad.
				shadowed = (trShad) ? scene->isShadowed(state, lightRay, sDepth, scol) : scene->isShadowed(state, lightRay);
				
				if(!shadowed && ls.pdf > 1e-6f)
				{
					if(trShad) ls.col *= scol;
					color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay);
					ls.col *= transmitCol;
					color_t surfCol = material->eval(state, sp, wo, lightRay.dir, BSDF_ALL);
					if( canIntersect)
					{
						float mPdf = material->pdf(state, sp, wo, lightRay.dir, BSDF_GLOSSY | BSDF_DIFFUSE | BSDF_DISPERSIVE | BSDF_REFLECT | BSDF_TRANSMIT);
						if(mPdf > 1e-6f)
						{
							float l2 = ls.pdf * ls.pdf;
							float m2 = mPdf * mPdf;
							float w = l2 / (l2 + m2);
							ccol += surfCol * ls.col * std::fabs(sp.N*lightRay.dir) * w / ls.pdf;
						}
						else
						{
							ccol += surfCol * ls.col * std::fabs(sp.N*lightRay.dir) / ls.pdf;
						}
					}
					else ccol += surfCol * ls.col * std::fabs(sp.N*lightRay.dir) / ls.pdf;
				}
			}
		}

		col += ccol * invNS;

		if(canIntersect) // sample from BSDF to complete MIS
		{
			color_t ccol2(0.f);

			hal2.setStart(offs-1);
			hal3.setStart(offs-1);

			for(int i=0; i<n; ++i)
			{
				ray_t bRay;
				bRay.tmin = MIN_RAYDIST; bRay.from = sp.P;

				float s1 = hal2.getNext();
				float s2 = hal3.getNext();

				sample_t s(s1, s2, BSDF_GLOSSY | BSDF_DIFFUSE | BSDF_DISPERSIVE | BSDF_REFLECT | BSDF_TRANSMIT);
				color_t surfCol = material->sample(state, sp, wo, bRay.dir, s);
				if( s.pdf>1e-6f && light->intersect(bRay, bRay.tmax, lcol, lightPdf) )
				{
					shadowed = (trShad) ? scene->isShadowed(state, bRay, sDepth, scol) : scene->isShadowed(state, bRay);
					if(!shadowed && lightPdf > 1e-6f)
					{
						if(trShad) lcol *= scol;
						color_t transmitCol = scene->volIntegrator->transmittance(state, lightRay);
						lcol *= transmitCol;
						float lPdf = 1.f/lightPdf;
						float l2 = lPdf * lPdf;
						float m2 = s.pdf * s.pdf;
						float w = m2 / (l2 + m2);
						CFLOAT cos2 = std::fabs(sp.N*bRay.dir);
						ccol2 += surfCol * lcol * cos2 * w / s.pdf;
					}
				}
			}
			col += ccol2 * invNS;
		}
	}
	return col;
}

bool mcIntegrator_t::createCausticMap()
{
	causticMap.clear();
	ray_t ray;
	std::vector<light_t *> causLights;
	
	for(unsigned int i=0;i<lights.size();++i)
	{
		if(lights[i]->shootsCausticP())
		{
			causLights.push_back(lights[i]);
		}
		
	}

	int numLights = causLights.size();
	progressBar_t *pb;
	if(intpb) pb = intpb;
	else pb = new ConsoleProgressBar_t(80);

	if(numLights > 0)
	{
		float lightNumPdf, lightPdf, s1, s2, s3, s4, s5, s6, s7, sL;
		float fNumLights = (float)numLights;
		float *energies = new float[numLights];
		for(int i=0;i<numLights;++i) energies[i] = causLights[i]->totalEnergy().energy();
		pdf1D_t *lightPowerD = new pdf1D_t(energies, numLights);
		
		Y_INFO << integratorName << ": Light(s) photon color testing for caustics map:" << yendl;
		color_t pcol(0.f);
		
		for(int i=0;i<numLights;++i)
		{
			pcol = causLights[i]->emitPhoton(.5, .5, .5, .5, ray, lightPdf);
			lightNumPdf = lightPowerD->func[i] * lightPowerD->invIntegral;
			pcol *= fNumLights*lightPdf/lightNumPdf; //remember that lightPdf is the inverse of the pdf, hence *=...
			Y_INFO << integratorName << ": Light [" << i+1 << "] Photon col:" << pcol << " | lnpdf: " << lightNumPdf << yendl;
		}
		
		delete[] energies;

		int pbStep;
		Y_INFO << integratorName << ": Building caustics photon map..." << yendl;
		pb->init(128);
		pbStep = std::max(1U, nCausPhotons / 128);
		pb->setTag("Building caustics photon map...");

		bool done=false;
		unsigned int curr=0;
		surfacePoint_t sp1, sp2;
		surfacePoint_t *hit=&sp1, *hit2=&sp2;
		renderState_t state;
		unsigned char userdata[USER_DATA_SIZE+7];
		state.userdata = (void *)( &userdata[7] - ( ((size_t)&userdata[7])&7 ) ); // pad userdata to 8 bytes
		while(!done)
		{
			state.chromatic = true;
			state.wavelength = RI_S(curr);
			s1 = RI_vdC(curr);
			s2 = scrHalton(2, curr);
			s3 = scrHalton(3, curr);
			s4 = scrHalton(4, curr);
			
			sL = float(curr) / float(nCausPhotons);
			
			int lightNum = lightPowerD->DSample(sL, &lightNumPdf);
			
			if(lightNum >= numLights)
			{
				Y_ERROR << integratorName << ": lightPDF sample error! " << sL << "/" << lightNum << yendl;
				delete lightPowerD;
				return false;
			}
			
			color_t pcol = causLights[lightNum]->emitPhoton(s1, s2, s3, s4, ray, lightPdf);
			ray.tmin = MIN_RAYDIST;
			ray.tmax = -1.0;
			pcol *= fNumLights * lightPdf / lightNumPdf; //remember that lightPdf is the inverse of th pdf, hence *=...
			if(pcol.isBlack())
			{
				++curr;
				done = (curr >= nCausPhotons);
				continue;
			}
			BSDF_t bsdfs = BSDF_NONE;
			int nBounces = 0;
			bool causticPhoton = false;
			bool directPhoton = true;
			const material_t *material = 0;
			const volumeHandler_t *vol = 0;
			
			while( scene->intersect(ray, *hit2) )
			{
				if(isnan(pcol.R) || isnan(pcol.G) || isnan(pcol.B))
				{
					Y_WARNING << integratorName << ": NaN (photon color)" << yendl;
					break;
				}
				color_t transm(1.f), vcol;
				// check for volumetric effects
				if(material)
				{
					if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(hit->Ng * ray.dir < 0)))
					{
						vol->transmittance(state, ray, vcol);
						transm = vcol;
					}

				}
				std::swap(hit, hit2);
				vector3d_t wi = -ray.dir, wo;
				material = hit->material;
				material->initBSDF(state, *hit, bsdfs);
				if(bsdfs & (BSDF_DIFFUSE | BSDF_GLOSSY))
				{
					//deposit caustic photon on surface
					if(causticPhoton)
					{
						photon_t np(wi, hit->P, pcol);
						causticMap.pushPhoton(np);
						causticMap.setNumPaths(curr);
					}
				}
				// need to break in the middle otherwise we scatter the photon and then discard it => redundant
				if(nBounces == causDepth) break;
				// scatter photon
				int d5 = 3*nBounces + 5;
				//int d6 = d5 + 1;

				s5 = scrHalton(d5, curr);
				s6 = scrHalton(d5+1, curr);
				s7 = scrHalton(d5+2, curr);

				pSample_t sample(s5, s6, s7, BSDF_ALL_SPECULAR | BSDF_GLOSSY | BSDF_FILTER | BSDF_DISPERSIVE, pcol, transm);
				bool scattered = material->scatterPhoton(state, *hit, wi, wo, sample);
				if(!scattered) break; //photon was absorped.
				pcol = sample.color;
				// hm...dispersive is not really a scattering qualifier like specular/glossy/diffuse or the special case filter...
				causticPhoton = ((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_DISPERSIVE)) && directPhoton) ||
								((sample.sampledFlags & (BSDF_GLOSSY | BSDF_SPECULAR | BSDF_FILTER | BSDF_DISPERSIVE)) && causticPhoton);
				// light through transparent materials can be calculated by direct lighting, so still consider them direct!
				directPhoton = (sample.sampledFlags & BSDF_FILTER) && directPhoton;
				// caustic-only calculation can be stopped if:
				if(!(causticPhoton || directPhoton)) break;

				if(state.chromatic && (sample.sampledFlags & BSDF_DISPERSIVE))
				{
					state.chromatic=false;
					color_t wl_col;
					wl2rgb(state.wavelength, wl_col);
					pcol *= wl_col;
				}
				ray.from = hit->P;
				ray.dir = wo;
				ray.tmin = MIN_RAYDIST;
				ray.tmax = -1.0;
				++nBounces;
			}
			++curr;
			if(curr % pbStep == 0) pb->update();
			done = (curr >= nCausPhotons);
		}
		pb->done();
		pb->setTag("Caustic photon map built.");
		Y_INFO << integratorName << ": Done." << yendl;
		Y_INFO << integratorName << ": Shot " << curr << " caustic photons from " << numLights <<" light(s)." << yendl;
		Y_INFO << integratorName << ": Stored caustic photons: " << causticMap.nPhotons() << yendl;
		
		delete lightPowerD;
		
		if(causticMap.nPhotons() > 0)
		{
			pb->setTag("Building caustic photons kd-tree...");
			causticMap.updateTree();
			Y_INFO << integratorName << ": Done." << yendl;
		}
		
		if(!intpb) delete pb;

	}
	else
	{
		Y_INFO << integratorName << ": No caustic source lights found, skiping caustic map building..." << yendl;
	}

	return true;
}

inline color_t mcIntegrator_t::estimateCausticPhotons(renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo) const
{
	if(!causticMap.ready()) return color_t(0.f);
	
	foundPhoton_t *gathered = new foundPhoton_t[nCausSearch];//(foundPhoton_t *)alloca(nCausSearch * sizeof(foundPhoton_t));
	int nGathered = 0;
	
	float gRadiusSquare = causRadius * causRadius;
	
	nGathered = causticMap.gather(sp.P, gathered, nCausSearch, gRadiusSquare);
	
	gRadiusSquare = 1.f / gRadiusSquare; 
	
	color_t sum(0.f);
	
	if(nGathered > 0)
	{
		const material_t *material = sp.material;
		color_t surfCol(0.f);
		float k = 0.f;
		const photon_t *photon;

		for(int i=0; i<nGathered; ++i)
		{
			photon = gathered[i].photon;
			surfCol = material->eval(state, sp, wo, photon->direction(), BSDF_ALL);
			k = kernel(gathered[i].distSquare, gRadiusSquare);
			sum += surfCol * k * photon->color();
		}
		sum *= 1.f / ( float(causticMap.nPaths()) );
	}
	
	delete [] gathered;
	
	return sum;
}

inline void mcIntegrator_t::recursiveRaytrace(renderState_t &state, diffRay_t &ray, BSDF_t bsdfs, surfacePoint_t &sp, vector3d_t &wo, color_t &col, float &alpha) const
{
	const material_t *material = sp.material;
	spDifferentials_t spDiff(sp, ray);

    state.raylevel++;
    
	if(state.raylevel <= rDepth)
	{
		Halton hal2(2);
		Halton hal3(3);

		// dispersive effects with recursive raytracing:
		if( (bsdfs & BSDF_DISPERSIVE) && state.chromatic)
		{
			state.includeLights = false; //debatable...
			int dsam = 8;
			int oldDivision = state.rayDivision;
			int oldOffset = state.rayOffset;
			float old_dc1 = state.dc1, old_dc2 = state.dc2;
			if(state.rayDivision > 1) dsam = std::max(1, dsam/oldDivision);
			state.rayDivision *= dsam;
			int branch = state.rayDivision*oldOffset;
			float d_1 = 1.f/(float)dsam;
			float ss1 = RI_S(state.pixelSample + state.samplingOffs);
			color_t dcol(0.f), vcol(1.f);
			vector3d_t wi;
			const volumeHandler_t *vol;
			diffRay_t refRay;
			
			for(int ns=0; ns<dsam; ++ns)
			{
				state.wavelength = (ns + ss1)*d_1;
				state.dc1 = scrHalton(2*state.raylevel+1, branch + state.samplingOffs);
				state.dc2 = scrHalton(2*state.raylevel+2, branch + state.samplingOffs);
				if(oldDivision > 1)  state.wavelength = addMod1(state.wavelength, old_dc1);
				state.rayOffset = branch;
				++branch;
				sample_t s(0.5f, 0.5f, BSDF_REFLECT|BSDF_TRANSMIT|BSDF_DISPERSIVE);
				color_t mcol = material->sample(state, sp, wo, wi, s);
				
				if(s.pdf > 1.0e-6f && (s.sampledFlags & BSDF_DISPERSIVE))
				{
					mcol *= std::fabs(wi*sp.N)/s.pdf;
					state.chromatic = false;
					color_t wl_col;
					wl2rgb(state.wavelength, wl_col);
					refRay = diffRay_t(sp.P, wi, MIN_RAYDIST);
					dcol += (color_t)integrate(state, refRay) * mcol * wl_col;
					state.chromatic = true;
				}
			}
			
			if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
			{
				vol->transmittance(state, refRay, vcol);
				dcol *= vcol;
			}
			
			col += dcol * d_1;
			state.rayDivision = oldDivision;
			state.rayOffset = oldOffset;
			state.dc1 = old_dc1;
			state.dc2 = old_dc2;
		}

		// glossy reflection with recursive raytracing:
		if( bsdfs & BSDF_GLOSSY )
		{
			state.includeLights = false;
			int gsam = 8;
			int oldDivision = state.rayDivision;
			int oldOffset = state.rayOffset;
			float old_dc1 = state.dc1, old_dc2 = state.dc2;
			if(state.rayDivision > 1) gsam = std::max(1, gsam/oldDivision);
			state.rayDivision *= gsam;
			int branch = state.rayDivision*oldOffset;
			unsigned int offs = gsam * state.pixelSample + state.samplingOffs;
			float d_1 = 1.f/(float)gsam;
			color_t gcol(0.f), vcol(1.f);
			vector3d_t wi;
			const volumeHandler_t *vol;
			diffRay_t refRay;

			hal2.setStart(offs);
			hal3.setStart(offs);

			for(int ns=0; ns<gsam; ++ns)
			{
				state.dc1 = scrHalton(2*state.raylevel+1, branch + state.samplingOffs);
				state.dc2 = scrHalton(2*state.raylevel+2, branch + state.samplingOffs);
				state.rayOffset = branch;
				++offs;
				++branch;

				float s1 = hal2.getNext();//RI_vdC(offs);
				float s2 = hal3.getNext();//scrHalton(2, offs);

				sample_t s(s1, s2, BSDF_REFLECT|BSDF_TRANSMIT|BSDF_FILTER|BSDF_GLOSSY);
				color_t mcol = material->sample(state, sp, wo, wi, s);

				if(s.pdf > 1.0e-6f && (s.sampledFlags & BSDF_GLOSSY))
				{
					mcol *= std::fabs(wi*sp.N)/s.pdf;
					refRay = diffRay_t(sp.P, wi, MIN_RAYDIST);
					gcol += (color_t)integrate(state, refRay) * mcol;
				}

				if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
				{
					if(vol->transmittance(state, refRay, vcol)) gcol *= vcol;
				}
			}

			col += gcol * d_1;
			state.rayDivision = oldDivision;
			state.rayOffset = oldOffset;
			state.dc1 = old_dc1;
			state.dc2 = old_dc2;
		}

		//...perfect specular reflection/refraction with recursive raytracing...
		if(bsdfs & (BSDF_SPECULAR | BSDF_FILTER) && state.raylevel < 20)
		{
			state.includeLights = true;
			bool reflect=false, refract=false;
			vector3d_t dir[2];
			color_t rcol[2], vcol;
			material->getSpecular(state, sp, wo, reflect, refract, dir, rcol);
			const volumeHandler_t *vol;
			if(reflect)
			{
				diffRay_t refRay(sp.P, dir[0], MIN_RAYDIST);
				spDiff.reflectedRay(ray, refRay);
				color_t integ = integrate(state, refRay);
				
				if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
				{
					if(vol->transmittance(state, refRay, vcol)) integ *= vcol;
				}
				
				col += integ * rcol[0];
			}
			if(refract)
			{
				diffRay_t refRay(sp.P, dir[1], MIN_RAYDIST);
				spDiff.refractedRay(ray, refRay, material->getMatIOR());
				colorA_t integ = integrate(state, refRay);
				
				if((bsdfs&BSDF_VOLUMETRIC) && (vol=material->getVolumeHandler(sp.Ng * refRay.dir < 0)))
				{
					if(vol->transmittance(state, refRay, vcol)) integ *= vcol;
				}
				
				col += (color_t)integ * rcol[1];
				alpha = integ.A;
			}
		}

	}
	--state.raylevel;
}

color_t mcIntegrator_t::sampleAmbientOcclusion(renderState_t &state, const surfacePoint_t &sp, const vector3d_t &wo) const
{
	color_t col(0.f), surfCol(0.f), scol(0.f);
	bool shadowed;
	const material_t *material = sp.material;
	ray_t lightRay;
	lightRay.from = sp.P;
	
	int n = aoSamples;
	if(state.rayDivision > 1) n = std::max(1, n / state.rayDivision);

	unsigned int offs = n * state.pixelSample + state.samplingOffs;
	
	Halton hal2(2);
	Halton hal3(3);

	hal2.setStart(offs-1);
	hal3.setStart(offs-1);

	for(int i = 0; i < n; ++i)
	{
		float s1 = hal2.getNext();
		float s2 = hal3.getNext();
		
		if(state.rayDivision > 1)
		{
			s1 = addMod1(s1, state.dc1);
			s2 = addMod1(s2, state.dc2);
		}
		
		lightRay.tmin = YAF_SHADOW_BIAS; // < better add some _smart_ self-bias value...this is still bad...
		lightRay.tmax = aoDist;		
		sample_t s(s1, s2, BSDF_GLOSSY | BSDF_DIFFUSE | BSDF_REFLECT );
		surfCol = material->sample(state, sp, wo, lightRay.dir, s);
		
		if(material->getFlags() & BSDF_EMIT)
		{
			col += material->emit(state, sp, wo) * s.pdf;
		}
		
		if(s.pdf > 1e-6f)
		{
			shadowed = (trShad) ? scene->isShadowed(state, lightRay, sDepth, scol) : scene->isShadowed(state, lightRay);
			
			if(!shadowed)
			{
				float cos = std::fabs(sp.N * lightRay.dir);
				if(trShad) col += aoCol * scol * surfCol * cos / s.pdf;
				else col += aoCol * surfCol * cos / s.pdf;
			}
		}
	}
	
	return col / (float)n;
}

void mcIntegrator_t::setICRecord(renderState_t &state, diffRay_t &ray, icRec_t *record) const {
	if (!ray.hasDifferentials)
		Y_INFO << "ERROR: ray from mcIntegrator_t::createNewICRecord() should have differentials" << std::endl;
	float oldRayLength[record->getM()];
	color_t oldRad[record->getM()];
	// we set the projected pixel area on the surface point
	record->setPixelArea(ray);
	ray_t sRay; // ray from hitpoint to hemisphere sample direction
	sRay.from = record->P;
	color_t innerRotValues;
	color_t innerTransValuesU;
	color_t innerTransValuesV;
	color_t radiance;
	unsigned int offs = 8 * state.pixelSample /*+ state.samplingOffs*/;
	record->stratHemi->randomize();
	for (int k=0; k<record->getN(); k++) {
		innerRotValues.black();
		innerTransValuesU.black();
		innerTransValuesV.black();
		for (int j=0; j<record->getM(); j++) {
			// Calculate each incoming radiance of hemisphere at point icRecord
			sRay.dir = record->getSampleHemisphere(j, k, offs);
			radiance = getRadiance(state, sRay);
			// note: oldRad[j] and oldRayLength[j] means L_j,k-1 and r_j,k-1 respectively
			//       oldRad[j-1] and oldRayLength[j-1] means L_j-1,k and r_j-1,k respectively
			if (k>0) {
				if (j>0) {
					// cos2(theta_j-)sin(theta_j-) * (L_j,k - L_j-1,k) / min(r_j,k , r_j-1,k)
					float cosThetaMin = record->stratHemi->getCosThetaMinus(j);
					innerTransValuesU +=
							( (cosThetaMin*cosThetaMin) * (radiance - oldRad[j-1]) )/
							(fmin(sRay.tmax, oldRayLength[j-1])) ;
					// cos(theta_j)[cos(theta_j-) - cos(theta_j+)] * (L_j,k - L_j,k-1) / [sin(theta_j,k) * min(r_j,k , r_j-1,k)]
					innerTransValuesV +=
							( record->stratHemi->getCosTheta(j) * (cosThetaMin - record->stratHemi->getCosThetaPlus(j)) *
								(radiance - oldRad[j]) ) /
							(record->stratHemi->getSinTheta(j) * fmin(sRay.tmax, oldRayLength[j]));
				}
			}
			record->irr += radiance;
			record->changeSampleRadius(sRay.tmax);
			// copy new rays and irradiance values over old ones
			oldRad[j] = radiance;
			oldRayLength[j] = sRay.tmax;
			innerRotValues -= record->stratHemi->getTanTheta(j) * radiance;
		}
		vector3d_t vk(record->stratHemi->getVk(k));
		record->rotGrad[0] += vk * innerRotValues.R;
		record->rotGrad[1] += vk * innerRotValues.G;
		record->rotGrad[2] += vk * innerRotValues.B;
		vector3d_t uk(record->stratHemi->getUk(k));
		vector3d_t vkm(record->stratHemi->getVkMinus(k));
		record->transGrad[0] +=
				( (innerTransValuesU.R * M_2PI / (float)record->getN() ) * uk ) +
				( innerTransValuesV.R * vkm );
		record->transGrad[1] +=
				( (innerTransValuesU.G * M_2PI / (float)record->getN() ) * uk ) +
				( innerTransValuesV.G * vkm );
		record->transGrad[2] +=
				( (innerTransValuesU.B * M_2PI / (float)record->getN() ) * uk ) +
				( innerTransValuesV.B * vkm );
	}
	float k = M_PI / ((float)record->getM() * (float)record->getN());
	record->irr = record->irr * k;
	for (int i=0; i<3; i++) {
		record->rotGrad[i] = changeBasis(
				record->rotGrad[i] * k,
				record->NU,
				record->NV,
				record->getNup());
		record->transGrad[i] = changeBasis(
				record->transGrad[i],
				record->NU,
				record->NV,
				record->getNup());
	}
	// HEURISTICS
	record->clampRbyGradient();
	if (ray.hasDifferentials)
		record->clampRbyScreenSpace();
	record->clampGradient();
}

void mcIntegrator_t::cleanup() {
	// do nothing, if IC implemented, may call the xml dump saving function
}

__END_YAFRAY
