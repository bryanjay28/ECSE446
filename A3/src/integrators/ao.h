/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Ambient occlusion integrator
 */
struct AOIntegrator : Integrator {

	// Use this in your switch statement to select the sampling type 
	ESamplingType m_samplingStrategy;

    explicit AOIntegrator(const Scene& scene) : Integrator(scene) { 
		m_samplingStrategy = scene.config.integratorSettings.ao.sampling_type;
	}

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);
		float L = 0.f;
		
		/*
		Use the m_sampling_type variable to set wi and the corresponding pdf 
		appropriately for sphere, hemisphere, or cosine sampling.

		You can use a switch statement or an if/else block.

		The m_sampling_type variable is an enum. The different values of the enum 
		can be accessed through:
		ESamplingType::ESpherical
		ESamplingType::EHemispherical
		ESamplingType::ECosineHemispherical
		*/
		
        // TODO(A3): Implement this
		// set variable for surface info
		SurfaceInteraction info;

		if (scene.bvh->intersect(ray, info)) {
			// sample a direction
			p2f randSample = sampler.next2D();

			// get the wi vector from the sampled direction vector using one of the warping strategies
			v3f wi(0.f);
			if (m_samplingStrategy == ESamplingType::ESpherical) {
				wi = Warp::squareToUniformSphere(randSample);
			}
			else if (m_samplingStrategy == ESamplingType::EHemispherical) {
				wi = Warp::squareToUniformHemisphere(randSample);
			}
			else if (m_samplingStrategy == ESamplingType::ECosineHemispherical) {
				wi = Warp::squareToCosineHemisphere(randSample);
			}
			wi = normalize(wi);

			// get the shadow ray, use the strategy from A2 
			// get the max distance of the shadow ray 
			float distance = scene.aabb.getBSphere().radius / 2;
			Ray shadow_ray = Ray(info.p, normalize(info.frameNs.toWorld(wi)), Epsilon, distance);
			
			// measure visibility in the AO integrand
			SurfaceInteraction shadowInfo;
			if (!scene.bvh->intersect(shadow_ray, shadowInfo)) {
				info.wi = wi;
				float cosTheta = dot(wi, normalize(info.frameNs.toLocal(info.frameNs.n)));

				// need to ensure it is visible
				if (Frame::cosTheta(wi) > 0 ) {
					// get the pdf using the associated samping strategy
					float pdf = 0.f;
					if (m_samplingStrategy == ESamplingType::ESpherical) {
						pdf = Warp::squareToUniformSpherePdf();
					}
					else if (m_samplingStrategy == ESamplingType::EHemispherical) {
						pdf = Warp::squareToUniformHemispherePdf(v3f(0));
					}
					else if (m_samplingStrategy == ESamplingType::ECosineHemispherical) {
						pdf = Warp::squareToCosineHemispherePdf(wi);
					}
					
					// Add up the total list base on the Ambient Oclussion formula
					// albedo is 1.0
					L = cosTheta / (M_PI * pdf);
				}
			}
			else L = 0.f;
		}
		Li = v3f(L);

        return Li;
    }
};

TR_NAMESPACE_END