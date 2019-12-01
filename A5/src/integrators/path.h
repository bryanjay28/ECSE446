/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Path tracer integrator
 */
struct PathTracerIntegrator : Integrator {
    explicit PathTracerIntegrator(const Scene& scene) : Integrator(scene) {
        m_isExplicit = scene.config.integratorSettings.pt.isExplicit;
        m_maxDepth = scene.config.integratorSettings.pt.maxDepth;
        m_rrDepth = scene.config.integratorSettings.pt.rrDepth;
        m_rrProb = scene.config.integratorSettings.pt.rrProb;
    }


    v3f renderImplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
        v3f Li(0.f), throughput(1.f);
        // TODO(A5): Implement this

		v3f emission;
		float pdf, totPdf = 1.f;

		// iterate until we reach the max number of bounces
		for (int depth = 0; depth <= m_maxDepth; depth++) {
			emission = getEmission(hit);
			
			// check if we hit the frontface of the emitter
			if (emission != v3f(0.f) && hit.wo.z > 0) {
				return emission * throughput / totPdf;
			}

			// calculate the bsdf to get the reflected vector and colour
			v3f bsdf = getBSDF(hit)->sample(hit, sampler, &pdf);

			if (pdf > 0) {
				throughput *= bsdf;
				totPdf *= pdf;
			}
			
			Ray rRay = Ray(hit.p + Epsilon, normalize(hit.frameNs.toWorld(hit.wi)));

			// if we don't find a surface then we return black as the ray never hits any object
			if (!scene.bvh->intersect(rRay, hit)) {
				return v3f(0.f);
			}
		}

        return Li;
    }

	v3f renderExplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& info) const {
		v3f Li(0.f), emission, L_bsdf(1.f);
		// TODO(A5): Implement this

		SurfaceInteraction hit = info;
		v3f throughput(1.f);

		//Checks if ray is on the emitter
		emission = getEmission(hit);

		if (emission != v3f(0.f)) {
			return emission;
		}
		
		// need to make an extra check for russian roulette which is when maxDepth is -1
		// is russian roulette we go until infinity and beyond
		for (int depth = 0; depth < m_maxDepth || m_maxDepth == -1; depth++) {
			// execute the surface area sampling like in direct A3

			// initialize values
			float emPdf, pdf;
			v3f ne, pos, wiW;

			// get the emitter 
			size_t shapeID = selectEmitter(sampler.next(), emPdf);
			const Emitter& em = getEmitterByID(shapeID);

			sampleEmitterPosition(sampler, em, ne, pos, pdf);

			// get wi and set wi to be in local coordinates
			wiW = normalize(pos - hit.p);
			hit.wi = hit.frameNs.toLocal(wiW);

			Ray sRay = Ray(hit.p, wiW);

			SurfaceInteraction sInfo;

			//Checks if there is no object intersection between the shading point and the light source
			if (scene.bvh->intersect(sRay, sInfo)) {
				emission = getEmission(sInfo);

				// verify if it lays on emitter and if shadow ray intersect with anything
				if (emission != v3f(0.f)) {
					// calcualte cosTheta to verify the visibility
					float cosTheta = dot(-wiW, ne);

					if (cosTheta > 0) {
						// find the jacobian and bsdf
						float jacobian = cosTheta / distance2(pos, hit.p);
						v3f brdf = getBSDF(hit)->eval(hit);

						Li += L_bsdf * jacobian * emission * brdf / (emPdf * pdf);
					}
				}
			}

			// execute the bsdf samplingf
			v3f throughput;

			// need to save the next hit point and set random emitter sample
			SurfaceInteraction nInfo;
			emission = v3f(1.f);

			// need to iterate until we hit the emitter
			while (emission != v3f(0.f)) {

				// calculate throughput and pdf like implicit
				throughput = getBSDF(hit)->sample(hit, sampler, &pdf);

				// create new ray from sample 
				Ray nRay = Ray(hit.p, normalize(hit.frameNs.toWorld(hit.wi)));

				// detect whether it intersects with another object if not return the Light
				if (!scene.bvh->intersect(nRay, nInfo)) {
					return Li;
				}

				// set new emission for while loop and set the next hit point
				emission = getEmission(nInfo);
			}

			hit = nInfo;
			// calculate the bsdf light
			if (pdf > 0) {
				L_bsdf *= throughput / pdf;
			}

			// do russian roulette path termination
			// we know m_rrDepth is 5 for 0-4 bounces in the TOML function so we will skip in explicit
			if (depth >= m_rrDepth) {
				// we use the probability to determine when the path will end
				if (sampler.next() > m_rrProb) {
					//	cout << "here" << depth << "\n";
					return Li;
				}
				L_bsdf /= m_rrProb;
			}

		}

		return Li;
	}

    v3f render(const Ray& ray, Sampler& sampler) const override {
        Ray r = ray;
        SurfaceInteraction hit;

        if (scene.bvh->intersect(r, hit)) {
            if (m_isExplicit)
                return this->renderExplicit(ray, sampler, hit);
            else
                return this->renderImplicit(ray, sampler, hit);
        }
        return v3f(0.0);
    }

    int m_maxDepth;     // Maximum number of bounces
    int m_rrDepth;      // When to start Russian roulette
    float m_rrProb;     // Russian roulette probability
    bool m_isExplicit;  // Implicit or explicit
};

TR_NAMESPACE_END
