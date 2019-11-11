/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once
#include <random>
#include "bsdfs/phong.h"

TR_NAMESPACE_BEGIN

/**
 * Reflection occlusion integrator
 */
struct ROIntegrator : Integrator {

    float m_exponent;

    explicit ROIntegrator(const Scene& scene) : Integrator(scene) {
        m_exponent = scene.config.integratorSettings.ro.exponent;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);
		float L = 0.f;

        // TODO(A3): Implement this
		// set variable for surface info
		SurfaceInteraction info;

		if (scene.bvh->intersect(ray, info)) {
			// sample a direction
			p2f randSample = sampler.next2D();

			// get the wi vector from the sampled direction vector using one of the warping strategies
			v3f wi = Warp::squareToPhongLobe(randSample, m_exponent);
			wi = normalize(wi);

			// get the pdf using the associated samping strategy using wi from phong frame
			float pdf = Warp::squareToPhongLobePdf(wi, m_exponent);

			// Get the reflection of wo so we can create a new fram around it
			v3f wr = reflect(info.wo);

			// Put vectors into the new frame
			Frame frame_R(wr);

			wi = frame_R.toWorld(wi);

			// calculate cosAlpha from dot product
			float cosAlpha = dot(normalize(wr), wi);	

			// get the shadow ray, use the strategy from A2 
			Ray shadow_ray = Ray(info.p, normalize(info.frameNs.toWorld(wi)));

			// need to ensure it is visible
			SurfaceInteraction shadowInfo;
			if (!scene.bvh->intersect(shadow_ray, shadowInfo)) {
				// get cosTheta the same way as before
				float cosTheta = dot(wi, normalize(info.frameNs.toLocal(info.frameNs.n)));
				
				// Add up the total list base on the reflective Oclussion formula given for Phongs
				L = (m_exponent + 2) * max(0.f, pow(cosAlpha, m_exponent)) * max(0.f, cosTheta) / (2 * M_PI * pdf);
			}
			else L = 0.f;
		}
		Li = v3f(L);

        return Li;
    }
};

TR_NAMESPACE_END