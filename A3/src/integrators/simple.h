/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Simple direct illumination integrator.
 */
struct SimpleIntegrator : Integrator {
    explicit SimpleIntegrator(const Scene& scene) : Integrator(scene) { }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);
		SurfaceInteraction info;

		// get the position and intensity of the light source
		v3f position = scene.getFirstLightPosition();
		v3f intensity = scene.getFirstLightIntensity();

		// Check the ray intersection similarly to A1 
		if (scene.bvh->intersect(ray, info)) {
			float distance = glm::distance(position, info.p);
			v3f light_pos = position - info.p;

			// get the vector pointing to the light (incident)
			info.wi = info.frameNs.toLocal(light_pos);

			// Create the shadow ray first between Epsilon (from ray definition) and distance to light
			Ray shadow_ray = Ray(info.p, normalize(light_pos), Epsilon, distance);
			
			if (scene.bvh->intersect(shadow_ray, info)) {
				return v3f(0.f);
			}

			// use built in function to get the BSDF
			v3f BSDF = getBSDF(info)->eval(info);

			Li = intensity * BSDF / pow(distance, 2);
		}
        // TODO(A2): Implement this

        return Li;
    }
};

TR_NAMESPACE_END