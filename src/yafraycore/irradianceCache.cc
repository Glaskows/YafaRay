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

#if HAVE_XML
#include <libxml/parser.h>
#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>
#define MY_ENCODING "ISO-8859-1"
#endif

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
	r = std::numeric_limits<float>::max();
}

icRec_t::icRec_t(int m, float kappa, const surfacePoint_t &sp):surfacePoint_t(sp), stratHemi(m), kappa(kappa) {
	r = std::numeric_limits<float>::max();
}

vector3d_t icRec_t::getSampleHemisphere(int j, int k) {
	return changeBasis( stratHemi.getDirection(j, k), NU, NV, Nup );
}

void icRec_t::changeSampleRadius(float newr) {
	// we use minimal distance radius (without clamping for now)
	if (newr < r) {
		r = newr;
	}
}

float icRec_t::getWeight(const icRec_t &record) const {
	float dot = Nup * record.getNup();
	// if record is pointing to the other side, better not to count his contribution
	if (dot<0.f) {
		return 0.f;
	}
	float epNor = fSqrt(1.f - dot) * NORMALIZATION_TERM;
	float epPos = (P - record.P).length() * 2.f / rClamp;
	float weight = 1.f - kappa * fmax(epPos, epNor);
	return weight;
}

bound_t icRec_t::getBound() const {
	return bound_t(P - rClamp, P + rClamp);
}

void icRec_t::setPixelArea(const diffRay_t &ray) {
	spDifferentials_t diff(*this, ray);
	pArea = fSqrt(diff.projectedPixelArea());
	rMin = 3.f * pArea;
	rMax = 20.f * pArea;
}

void icRec_t::setNup(const vector3d_t &wo) {
	Nup = FACE_FORWARD(Ng, N, wo);
}

bool icRec_t::inFront(const icRec_t &record) const {
	float di = (P - record.P) * ((Nup + record.getNup() )/2.0f);
	if (di < -0.01f) // small negative value, ¿it works?
		return true;
	return false;
}

void icRec_t::clampRbyGradient() {
	rClamp = std::min(r, std::min(irr.R/transGrad[0].length(),
								  std::min( irr.G/transGrad[1].length(), irr.B/transGrad[2].length())) );
}

void icRec_t::clampRbyScreenSpace() {
	rClamp = std::min(std::max(rClamp, rMin), rMax);
}

void icRec_t::clampRbyNeighbor() {
	rClamp = std::min(rClamp, rNeighbor);
}

void icRec_t::clampGradient() {
	for (int i=0; i<3; i++) {
		transGrad[i] =
				transGrad[i] * std::min(1.f, r/rMin);
	}
}

void icRec_t::setRNeighbor(float r) {
	rNeighbor = r;
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
	if (!record.inFront(sample)) {
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
	} else {
	//	Y_INFO << "In front!" << std::endl;
	}
	return true; // when could it be false? example?
}

void icTree_t::recursiveFindNear(octNode_t<icRec_t> *node, const bound_t &nodeBound,
								 const icRec_t &record, std::vector<icRec_t *> &nearRecs,
								 float &minR) {
	for (unsigned int i = 0; i < node->data.size(); ++i) {
		// pass the "in front" test
		if (!record.inFront(node->data[i])) {
			float distance = (record.P - (node->data[i]).P).length();
			// if both radius overlaps
			if ( distance <= (record.r + (node->data[i]).r) ) {
				float rSum = (node->data[i]).r + distance;
				// checks for triangule inequality
				if ( rSum < record.r ){
					minR = std::min( minR, rSum );
					// add pointer to record to nearRecs
					nearRecs.push_back(&(node->data[i]));
				}
			}
		}
	}
	bound_t dataBound(record.P-record.r, record.P+record.r);
	// check on all the childrens that the records radius overlap
	point3d_t center = nodeBound.center();
	// Determine which children the item overlaps
	bool over[8];
	over[1] = over[3] = over[5] = over[7] = (dataBound.a.x <= center.x);
	over[0] = over[2] = over[4] = over[6] = (dataBound.g.x  > center.x);
	if(dataBound.a.y > center.y)  over[2] = over[3] = over[6] = over[7] = false;
	if(dataBound.g.y <= center.y) over[0] = over[1] = over[4] = over[5] = false;
	if(dataBound.a.z > center.z)  over[4] = over[5] = over[6] = over[7] = false;
	if(dataBound.g.z <= center.z) over[0] = over[1] = over[2] = over[3] = false;
	for (int child = 0; child < 8; ++child)
	{
		// don't do anything if radius not overlap or child node doent exist
		if (!over[child] || !node->children[child]) continue;
		// compute childbound and keep searching in child
		bound_t childBound;
		childBound.a.x = (child & 1) ? nodeBound.a.x : center.x;
		childBound.g.x = (child & 1) ? center.x : nodeBound.g.x;
		childBound.a.y = (child & 2) ? nodeBound.a.y : center.y;
		childBound.g.y = (child & 2) ? center.y : nodeBound.g.y;
		childBound.a.z = (child & 4) ? nodeBound.a.z : center.z;
		childBound.g.z = (child & 4) ? center.z : nodeBound.g.z;
		recursiveFindNear(node->children[child], childBound, record, nearRecs, minR);
	}
}

void icTree_t::neighborClamp(icRec_t &record) {
	// perform lock in tree, so we can't add other records at the same time
	lock.readLock();
	// if record is outside the scene don't do anything
	if (!treeBound.includes(record.P)) return;
	// create the vector to store all the neighbors
	std::vector<icRec_t *> nearRecs;
	// set neighbor clamped radius equal to the original distance to surfaces
	float minR = record.r;
	// search for all the near records
	recursiveFindNear(&root, treeBound, record, nearRecs, minR);
	// perform neighbor clamp on record
	record.setRNeighbor(minR);
	record.clampRbyNeighbor();
	// perform neighbor clamp to neighbors
	for (unsigned int i=0; i<nearRecs.size(); i++) {
		nearRecs[i]->setRNeighbor(minR + (nearRecs[i]->P - record.P).length());
		nearRecs[i]->clampRbyNeighbor();
	}
	// TODO: neighbor clamp algorith
	lock.unlock();
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

void icTree_t::saveToXml(const std::string &fileName) {
	int rc;
	xmlTextWriterPtr writer;
	xmlChar *tmp;
	Y_INFO << "Start dump of IC Tree to xml file: " << fileName << std::endl;
	// Create a new XmlWriter for uri
	writer = xmlNewTextWriterFilename(fileName.c_str(), 0);
	if (writer == NULL) {
		Y_INFO << "testXmlwriterFilename: Error creating the xml writer" << std::endl;
		return;
	}
	// Start the document
	rc = xmlTextWriterStartDocument(writer, NULL, MY_ENCODING, NULL);
	if (rc < 0) {
		Y_INFO << "testXmlwriterFilename: Error at xmlTextWriterStartDocument" << std::endl;
		return;
	}
	// Create root element named ICtree
	rc = xmlTextWriterStartElement(writer, BAD_CAST "ICtree");
	if (rc < 0) {
		Y_INFO << "testXmlwriterFilename: Error at xmlTextWriterStartElement\n" << std::endl;
		return;
	}

	// RECURSIVE? saveNodeToXml(writer, & root);

	octNode_t<icRec_t> *nodes[maxDepth+1];
	int child[maxDepth+1];

	int level = 0;
	nodes[0] = &root;
	child[0] = 0;

	do {
		// process data
		Y_INFO << "level: " << level << " - Sibling: " << child[level] << " - N° nodos: " << nodes[level]->data.size() << std::endl;

		// go down one level down
		level++;

		nodes[level] = nodes[level-1]->children[child[level-1]]; //nodes[level-1]->children[0];
		child[level] = 0;

		// until we find a valid node
		while (!nodes[level]) {
			// try to go to next ancestor brother
			do {
				level--;
			} while (child[level]==7 && level>=0);
			if (level >= 0) {
				child[level]++;
				level++;
				nodes[level] = nodes[level-1]->children[child[level-1]];
				child[level] = 0;
			} else {
				break; // get out of valid node check
			}
		}
	} while (level >= 0 );

	// Close the element named ICtree
	rc = xmlTextWriterEndElement(writer);
	if (rc < 0) {
		Y_INFO << "testXmlwriterFilename: Error at xmlTextWriterEndElement\n" << std::endl;
		return;
	}

	xmlFreeTextWriter(writer);
	Y_INFO << "End xml dump" << std::endl;
}


__END_YAFRAY
