/*
	This file is part of TinyRender, an educative rendering system.
	Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
	Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <tiny_obj_loader.h>
#define RAY_EPS_CV 1e-5 // Use when setting min and max dist for ray in control variates code

TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator for polygonal light sources
 * Follows Arvo '94.
 */
struct PolygonalIntegrator : Integrator {

	float m_alpha;             // Control variates "strength"
	size_t m_visSamples;       // # of samples to estimate h - alpha*g
	bool m_traceShadows;       // Trace shadows or not
	EPolygonalMethod m_method; // Method to use (Arvo, or control variates)

	std::vector<std::vector<v3f>> m_triangles; // Data structure to store triangles

	explicit PolygonalIntegrator(const Scene& scene) : Integrator(scene) {
		m_alpha = scene.config.integratorSettings.poly.alpha;
		m_visSamples = scene.config.integratorSettings.poly.visSamples;
		m_traceShadows = scene.config.integratorSettings.poly.traceShadows;
		m_method = scene.config.integratorSettings.poly.method;

		/**
		 * 1) Get # of triangles on emitter
		 * 2) Store vertices in m_triangles
		 */
		 // TODO(A4): Implement this

		size_t shapeID = scene.getFirstLight();
		auto shape = scene.worldData.shapes[shapeID];

		int mesh_len = shape.mesh.indices.size();

		for (int i = 0; i < shape.mesh.num_face_vertices.size(); i++) {
			int id = i * 3;

			vector<v3f> tri(3);
			tri[0] = scene.getObjectVertexPosition(shapeID, id);
			tri[1] = scene.getObjectVertexPosition(shapeID, id + 1);
			tri[2] = scene.getObjectVertexPosition(shapeID, id + 2);

			m_triangles.push_back(tri);
		}
	}

	/// Reflect
	inline v3f reflect(const v3f& d) const {
		return v3f(-d.x, -d.y, d.z);
	}

	/**
	 * === PHONG BONUS ONLY ===
	 * Compute the following integral:
	 *    T(a, b, n, x) = \int_0^x [a \cos(\theta) + b \sin(\theta)]ˆn d\theta
	 * Uses a recurrent relation (see Snyder's note, 1996)
	 *
	 * Series function:
	 *    T_sum(a, b, n, x) = \sum_{i=0}ˆ{(n-1)/2} T(a, b, 2i+1, x)
	 * assuming n is _odd_
	 */
	float cosineSinePowerIntegralSum(float a, float b, int exp, float theta) const {
		if (exp % 2 == 0) exp += 1; // Make exponent odd
		float Tsum = 0.f;

		// Implementing this function may be useful if you attempt the bonus

		// TODO(A4): Implement this

		return Tsum;
	}

	/**
	 * Compute edge (v1--v2) contribution
	 * The exp term is only needed if you attempt the bonus, otherwise, you can ignore it
	 */
	float getEdgeContrib(const v3f& v1, const v3f& v2, const SurfaceInteraction& i, int exp = 0) const {
		float contrib = 0.f;
		// TODO(A4): Implement this

		v3f temp1 = normalize(v1 - i.p);
		v3f temp2 = normalize(v2 - i.p);
		
		v3f gamma_i =  normalize(cross(temp1, temp2));
		float theta = acos(dot(temp1, temp2));
		
		contrib = theta * dot(-gamma_i, i.frameNs.n);

		return contrib;
	}


	/// Direct illumination using Arvo '94 analytic solution for polygonal lights
	v3f renderAnalytic(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		// TODO(A4): Implement this
		SurfaceInteraction info;

		if (scene.bvh->intersect(ray, info)) {
			v3f emission = getEmission(info);

			// check the emission radiance
			if (emission != v3f(0.f)) {
				return emission;
			}

			v3f E_poly(0.f);
			float sum = 0.f;

			// calculate the sum of all the edge contributions of all the triangles
			for (int i = 0; i < m_triangles.size(); i++) {
				sum += getEdgeContrib(m_triangles[i][0], m_triangles[i][1], info);
				sum += getEdgeContrib(m_triangles[i][1], m_triangles[i][2], info);
				sum += getEdgeContrib(m_triangles[i][2], m_triangles[i][0], info);
			}

			size_t shapeID = scene.getFirstLight();
			const Emitter& em = getEmitterByID(getEmitterIDByShapeID(shapeID));

			// get E poly based on the eqn from the handout
			E_poly = em.getPower() * sum / (2.f * em.area * M_PI);

			info.wi = v3f(0, 0, 1.f);
			v3f brdf = getBSDF(info)->eval(info);

			Lr = brdf * E_poly;
		}

		return Lr;
	}

	/**
	 * Stand-alone estimator for h - alpha*g (with primary ray)
	 * Trace a primary ray, check for emitter hit, and then call `estimateVisDiff()`
	 * Used by polygonal render pass
	 */
	v3f estimateVisDiffRealTime(const Ray& ray, Sampler& sampler, const Emitter& em) {
		v3f D(0.f);

		SurfaceInteraction hit;
		if (!scene.bvh->intersect(ray, hit)) return D;

		const BSDF* bsdf = getBSDF(hit);
		if (bsdf->isEmissive()) return D;

		hit.wi = v3f(0, 0, 1); // Trick to get 1/pi * albedo without cosine term
		D = estimateVisDiff(sampler, hit, em);

		return D;
	}

	/// Stand-alone estimator for h - alpha*g (without primary ray)
	/// Use RAY_EPS_CV when setting min and max dist for shadow ray
	v3f estimateVisDiff(Sampler& sampler, SurfaceInteraction& i, const Emitter& em) const {
		v3f sum(0.f);
		// TODO(A4): Implement this

		v3f g(0.f), h(0.f);

		v3f pe;    // Point on emitter
		v3f ne;    // Surface normal at point
		float pdf; // PDF of choosing point;

		sampleEmitterPosition(sampler, em, ne, pe, pdf); // Sample mesh uniformly

		v3f emitterCenter = scene.getShapeCenter(em.shapeID);
		float emPdf = 1.f / scene.emitters.size();
		v3f pShading = i.p;

		v3f wiW = normalize(pe - pShading);
		i.wi = i.frameNs.toLocal(wiW);

		SurfaceInteraction sInfo; 
		Ray shadowRay = Ray(i.p, wiW, RAY_EPS_CV);

		if (scene.bvh->intersect(shadowRay, sInfo)) {
			v3f newEmission = getEmission(sInfo);
			if (newEmission != v3f(0.f)) {
				v3f visibilityVector = emitterCenter - pShading;

				float visibilityPt = dot(normalize(visibilityVector), -ne);

				if (visibilityPt >= 0) {
					// get the normal of wi in y direction corresponding to pdf
					float wi_normal_y = dot(wiW, -ne);

					// find the new pdf based on omega
					float pdf_omega = pdf * pow(distance(pe, pShading), 2) / wi_normal_y;

					// get the brdf
					v3f brdf = getBSDF(i)->eval(i);

					if (pdf_omega > 0) {
						// add up the Lr light then take hte avg later
						h += brdf * newEmission / (pdf_omega * emPdf);
					}
				}
			}
		}

		v3f visibilityVector = emitterCenter - pShading;
		float visibilityPt = dot(normalize(visibilityVector), -ne);

		if (visibilityPt >= 0) {
			// get the normal of wi in y direction corresponding to pdf
			float wi_normal_y = dot(wiW, -ne);

			// find the new pdf based on omega
			float pdf_omega = pdf * pow(distance(pe, pShading), 2) / wi_normal_y;

			// get the brdf
			v3f brdf = getBSDF(i)->eval(i);

			if (pdf_omega > 0) {
				// add up the Lr light then take hte avg later
				g += brdf * em.getRadiance() / (pdf_omega * emPdf);
			}
		}

		sum += h - m_alpha * g;
		

		return sum;
	}

	/// Control variates using Arvo '94 for direct illumination; ray trace shadows

	v3f renderControlVariates(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		// TODO(A4): Implement this

		SurfaceInteraction info;

		size_t shapeID = scene.getFirstLight();
		const Emitter& em = getEmitterByID(getEmitterIDByShapeID(shapeID));

		if (scene.bvh->intersect(ray, info)) {
			v3f emission = getEmission(info);

			if (emission != v3f(0.f)) {
				return emission;
			}

			// compute the render analytical here
			v3f E_poly(0.f);
			float total = 0.f;

			// calculate the sum of all the edge contributions of all the triangles
			for (int j = 0; j < m_triangles.size(); j++) {
				total += getEdgeContrib(m_triangles[j][0], m_triangles[j][1], info);
				total += getEdgeContrib(m_triangles[j][1], m_triangles[j][2], info);
				total += getEdgeContrib(m_triangles[j][2], m_triangles[j][0], info);
			}

			// get E poly based on the eqn from the handout
			E_poly = em.getPower() * total / (2.f * em.area * M_PI);

			info.wi = v3f(0, 0, 1.f);
			v3f brdf = getBSDF(info)->eval(info);

			Lr += brdf * E_poly;

			v3f controlVariates(0.f);

			for (int i = 0; i < m_visSamples; i++) {
				controlVariates += estimateVisDiff(sampler, info, em);
			}

			Lr += controlVariates / m_visSamples;
		}

		Lr = clampBelow(Lr, 0);

		return Lr;
	}

	/// Direct illumination using surface area sampling
	v3f renderArea(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		// TODO(A4): Implement this

		SurfaceInteraction info;

		if (scene.bvh->intersect(ray, info)) {
			// need to verify if intersection lies on surface of the emitter
			v3f emission = getEmission(info);

			if (emission != v3f(0.f)) {
				return emission;
			}

			v3f pe;    // Point on emitter
			v3f ne;    // Surface normal at point
			float pdf; // PDF of choosing point

			size_t shapeID = scene.getFirstLight();
			const Emitter& em = getEmitterByID(getEmitterIDByShapeID(shapeID));

			sampleEmitterPosition(sampler, em, ne, pe, pdf); // Sample mesh uniformly
				
			v3f emitterCenter = scene.getShapeCenter(em.shapeID);
			float emPdf = 1.f / scene.emitters.size();
			v3f pShading = info.p;

			// must check the visibility point
			v3f wiW = normalize(pe - pShading);
			info.wi = info.frameNs.toLocal(wiW);

			if (m_traceShadows) {
				SurfaceInteraction sInfo;
				Ray shadowRay = Ray(info.p, normalize(wiW));

				if (scene.bvh->intersect(shadowRay, sInfo)) {
					// check where the intersection lies again
					v3f newEmission = getEmission(sInfo);
					if (newEmission != v3f(0.f)) {
						v3f visibilityVector = emitterCenter - pShading;

						float visibilityPt = dot(normalize(visibilityVector), -ne);

						if (visibilityPt >= 0) {
							// get the normal of wi in y direction corresponding to pdf
							float wi_normal_y = dot(wiW, -ne);

							// find the new pdf based on omega
							float pdf_omega = pdf * pow(distance(pe, pShading), 2) / wi_normal_y;

							// get the brdf
							v3f brdf = getBSDF(info)->eval(info);

							if (pdf_omega > 0) {
								// add up the Lr light then take hte avg later
								Lr += brdf * newEmission / (pdf_omega * emPdf);
							}
						}
					}
				}					
			}
			else {
				v3f visibilityVector = emitterCenter - pShading;
				float visibilityPt = dot(normalize(visibilityVector), -ne);

				if (visibilityPt >= 0) {
					// get the normal of wi in y direction corresponding to pdf
					float wi_normal_y = dot(wiW, -ne);

					// find the new pdf based on omega
					float pdf_omega = pdf * pow(distance(pe, pShading), 2) / wi_normal_y;

					// get the brdf
					v3f brdf = getBSDF(info)->eval(info);

					if (pdf_omega > 0) {
						// add up the Lr light then take hte avg later
						Lr += brdf * em.getRadiance() / (pdf_omega * emPdf);
					}
				}
			}
		}

		return Lr;
	}

	/// Branch to corresponding method
	v3f render(const Ray& ray, Sampler& sampler) const override {
		switch (m_method) {
		case EPolygonalMethod::ESurfaceArea:
			return PolygonalIntegrator::renderArea(ray, sampler);
			break;
		case EPolygonalMethod::EControlVariates:
			return PolygonalIntegrator::renderControlVariates(ray, sampler);
			break;
		default:
			return PolygonalIntegrator::renderAnalytic(ray, sampler);
			break;
		}
	}

};

TR_NAMESPACE_END