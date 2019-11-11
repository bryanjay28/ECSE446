/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator with MIS
 */
struct DirectIntegrator : Integrator {
    explicit DirectIntegrator(const Scene& scene) : Integrator(scene) {
        m_emitterSamples = scene.config.integratorSettings.di.emitterSamples;
        m_bsdfSamples = scene.config.integratorSettings.di.bsdfSamples;
        m_samplingStrategy = scene.config.integratorSettings.di.samplingStrategy;
    }

    static inline float balanceHeuristic(float nf, float fPdf, float ng, float gPdf) {
        float f = nf * fPdf, g = ng * gPdf;
        return f / (f + g);
    }

    void sampleSphereByCosineHemisphere(const p2f& sample,
                                        const v3f& n,
                                        const p3f& pShading,
                                        const v3f& emitterCenter,
                                        float emitterRadius,
                                        v3f& wiW,
                                        float& pdf) const {
		v3f wi = Warp::squareToCosineHemisphere(sample);
		pdf = Warp::squareToCosineHemispherePdf(wi);

		Frame frameNs(n);
		wiW = frameNs.toWorld(wi);

        // TODO(A3): Implement this
    }

    void sampleSphereByArea(const p2f& sample,
                            const p3f& pShading,
                            const v3f& emitterCenter,
                            float emitterRadius,
                            v3f& pos,
                            v3f& ne,
                            v3f& wiW,
                            float& pdf) const {
        // TODO(A3): Implement this
		// pdf is the inverse area of a sphere
		pdf = 1 / (4 * M_PI * pow(emitterRadius, 2));

		v3f square_wi = Warp::squareToUniformSphere(sample);

		// get the position of the light emitter, needed to scale the output with radius
		pos = emitterCenter + emitterRadius * square_wi;
		ne = normalize(emitterRadius * square_wi);

		// must check the visibility point
		wiW = normalize(pos - pShading);

    }

    void sampleSphereBySolidAngle(const p2f& sample,
                                  const p3f& pShading,
                                  const v3f& emitterCenter,
                                  float emitterRadius,
                                  v3f& wiW,
                                  float& pdf) const {
        // TODO(A3): Implement this
		// calculate the costheta max which is jus the magnitude
		float dist = distance(emitterCenter, pShading);
		float cosThetaMax = dist / sqrt(pow(dist, 2) + pow(emitterRadius, 2));

		pdf = Warp::squareToUniformConePdf(cosThetaMax);
		v3f wi = Warp::squareToUniformCone(sample, cosThetaMax);

		// get the vector wc that we will will put our wi in its frame
		v3f wc = normalize(emitterCenter - pShading);
		Frame frame_wc(wc);

		wiW = normalize(frame_wc.toWorld(wi));
    }

    v3f renderArea(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        // TODO(A3): Implement this

		SurfaceInteraction info;

		if (scene.bvh->intersect(ray, info)) {
			// need to verify if intersection lies on surface of the emitter
			v3f emission = getEmission(info);

			if (emission == v3f(0)) {
				// perform the estimator around a sphere
				// set the parameters from the random emitter sample
				float emPdf, pdf;
				v3f wiW, pos, ne;
				v3f pShading = info.p;

				// iterate through based on the emitter samples
				for (int i = 0; i < m_emitterSamples; i++) {
					// get random sample
					v2f sample = sampler.next2D();

					size_t id = selectEmitter(sampler.next(), emPdf);
					const Emitter& em = getEmitterByID(id);
					float emitterRadius = scene.getShapeRadius(em.shapeID);
					v3f emitterCenter = scene.getShapeCenter(em.shapeID);

					// get the sample for wi and the pdf
					sampleSphereByArea(sample, pShading, emitterCenter, emitterRadius, pos, ne, wiW, pdf);
					info.wi = info.frameNs.toLocal(wiW);

					// trace the shadow ray and verify intersection
					SurfaceInteraction sInfo;
					Ray shadowRay = Ray(info.p, normalize(wiW));

					if (scene.bvh->intersect(shadowRay, sInfo)) {
						// check where the intersection lies again
						v3f newEmission = getEmission(sInfo);

						if (newEmission != v3f(0)) {
							// verify visibility of the light to see if its blocked
							v3f visibilityVector = emitterCenter - pShading;
							float visibilityPt = dot(normalize(visibilityVector), ne);

							if (visibilityPt >= 0) {
								// get the normal of wi in y direction corresponding to pdf
								float wi_normal_y = dot(wiW, ne);
								
								// fing the new pdf based on omega
								float pdf_omega = pdf * pow(distance(pos, pShading), 2) / wi_normal_y;

								// get the brdf
								v3f brdf = getBSDF(info)->eval(info);

								if (pdf_omega > 0) {
									// add up the Lr light then take hte avg later
									Lr += brdf * newEmission / (pdf_omega * emPdf);
								}
								else Lr += 0;

							}
						}
					}
				}
			} // return the emission since intersection is on the surface
			else return emission;
		}

		// take the average here before we return based on the number of emission samples 
		return Lr / m_emitterSamples;
    }

    v3f renderCosineHemisphere(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

		SurfaceInteraction info;
		
		// check for ray intersection
		if (scene.bvh->intersect(ray, info)) {
			// need to verify if intersection lies on surface of the emitter
			v3f emission = getEmission(info);

			if (emission == v3f(0)) {
				// perform monte carlo estimator
				float pdf, emitterRadius = 0.f;
				v3f pShading, emitterCenter, wiW;

				// iterate through based on the emitter samples
				for (int i = 0; i < m_bsdfSamples; i++) {
					// get random sample
					v2f sample = sampler.next2D();

					// get the sample for wi and the pdf
					sampleSphereByCosineHemisphere(sample, info.frameNs.n, pShading, emitterCenter, emitterRadius, wiW, pdf);
					info.wi = info.frameNs.toLocal(wiW);

					// trace the shadow ray and verify intersection
					SurfaceInteraction sInfo;
					Ray shadowRay = Ray(info.p, normalize(wiW));

					if (scene.bvh->intersect(shadowRay, sInfo)) {
						// check where the intersection lies again
						v3f newEmission = getEmission(sInfo);
						if (newEmission != v3f(0)) {
							// get the brdf
							v3f brdf = getBSDF(info)->eval(info);

							// add up the Lr light then take hte avg later
							if (pdf <= 0) {
								Lr += 0;
							}
							else Lr += brdf * newEmission / pdf;
						}
					}
				}
			} // return the emission since intersection is on the surface
			else return emission;
		}

        // TODO(A3): Implement this
		// take the average here before we return based on the number of emission samples 
		return Lr / m_bsdfSamples;
    }

    v3f renderBSDF(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

		SurfaceInteraction info;

		// check for ray intersection
		if (scene.bvh->intersect(ray, info)) {
			// need to verify if intersection lies on surface of the emitter
			v3f emission = getEmission(info);

			if (emission == v3f(0.f)) {
				// perform monte carlo estimator
				float pdf = 0.f;

				// iterate through based on the emitter samples
				for (int i = 0; i < m_bsdfSamples; i++) {
					// trace the shadow ray and verify intersection
					v3f brdf = getBSDF(info)->sample(info, sampler, &pdf);
					
					SurfaceInteraction sInfo;
					Ray shadowRay = Ray(info.p, normalize(info.frameNs.toWorld(info.wi)));

					if (scene.bvh->intersect(shadowRay, sInfo)) {
						// check where the intersection lies again
						emission = getEmission(sInfo);

						if (emission != v3f(0.f)) {
							// add up the Lr light then take hte avg later
							Lr += brdf * emission;
						}
					}
				}
			} // return the emission since intersection is on the surface
			else return emission;
		}

        // TODO(A3): Implement this
		// take avg
		return Lr / m_bsdfSamples;
    }

    v3f renderSolidAngle(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

		// perform similar method to before

		// set surface interaction 
		SurfaceInteraction info;

		if (scene.bvh->intersect(ray, info)) {
			// need to verify if intersection lies on surface of the emitter
			v3f emission = getEmission(info);

			if (emission == v3f(0)) {
				// perform the estimator around a sphere
				// set the parameters from the random emitter sample
				float emPdf, pdf;
				v3f wiW, pShading = info.p;

				// iterate through based on the emitter samples
				for (int i = 0; i < m_emitterSamples; i++) {
					// get random sample
					v2f sample = sampler.next2D();

					size_t id = selectEmitter(sampler.next(), emPdf);
					const Emitter& em = getEmitterByID(id);
					float emitterRadius = scene.getShapeRadius(em.shapeID);
					v3f emitterCenter = scene.getShapeCenter(em.shapeID);

					// get the sample for wi and the pdf
					sampleSphereBySolidAngle(sample, pShading, emitterCenter, emitterRadius, wiW, pdf);
					info.wi = info.frameNs.toLocal(wiW);

					// trace the shadow ray and verify intersection
					SurfaceInteraction sInfo;
					Ray shadowRay = Ray(info.p, normalize(wiW));

					if (scene.bvh->intersect(shadowRay, sInfo)) {
						// check where the intersection lies again
						v3f newEmission = getEmission(sInfo);

						if (newEmission != v3f(0)) {
							// get the brdf and sum up Lr
							v3f brdf = getBSDF(info)->eval(info);

							// ensure pdf is positive
							if (pdf > 0) {
								// add up the Lr light then take hte avg later
								Lr += brdf * newEmission / (pdf * emPdf);
							}
							else Lr += 0;
						}
					}
					else i--;
				}
			} // return the emission since intersection is on the surface
			else return emission;
		}

		// take the average over the emitter samples 
		return Lr / m_emitterSamples;
    }

    v3f renderMIS(const Ray& ray, Sampler& sampler) const {

        v3f Lr(0.f);

        // TODO(A4): Implement this

        return Lr;
    }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        if (m_samplingStrategy == ESamplingStrategy::EMIS)
            return this->renderMIS(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::EArea)
            return this->renderArea(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::ESolidAngle)
            return this->renderSolidAngle(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::ECosineHemisphere)
            return this->renderCosineHemisphere(ray, sampler);
        else
            return this->renderBSDF(ray, sampler);
    }

    size_t m_emitterSamples;     // Number of emitter samples
    size_t m_bsdfSamples;        // Number of BSDF samples
    ESamplingStrategy m_samplingStrategy;   // Sampling strategy to use
};

TR_NAMESPACE_END