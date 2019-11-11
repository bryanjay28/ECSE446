/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Perfectly diffuse, Lambertian reflectance model
 */
struct DiffuseBSDF : BSDF {
    std::unique_ptr<Texture < v3f>> albedo;

    DiffuseBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        const tinyobj::material_t& mat = scene.materials[matID];

        if (mat.diffuse_texname.empty())
            albedo = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
        else
            albedo = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

        components.push_back(EDiffuseReflection);

        combinedType = 0;
        for (size_t i = 0; i < components.size(); ++i)
            combinedType |= components[i];
    }

    v3f eval(const SurfaceInteraction& i) const override {
        v3f val(0.f);
		v3f black_colour(0.f);

		// need to normalize both vectors 
		v3f normal_wo = normalize(i.wo);
		v3f normal_wi = normalize(i.wi);

        // TODO(A2): Implement this

		v3f albedo_rho = albedo->eval(worldData, i);

		// Perform to check to see if rays are front facing to the surface.
		if (!(Frame::cosTheta(normal_wo) > 0 && Frame::cosTheta(normal_wi) > 0)) {
			return black_colour;
		}
		
		// return the albedo divided by π and multiplied by the cosine foreshortening factor 
		val = Frame::cosTheta(normal_wi) * albedo_rho / M_PI;
        return val;
    }

    float pdf(const SurfaceInteraction& i) const override {
        float pdf = 0.f;

		pdf = Warp::squareToCosineHemispherePdf(i.wi);

        // TODO(A3): Implement this
        return pdf;
    }

    v3f sample(SurfaceInteraction& i, Sampler& sampler, float* _pdf) const override {
        v3f val(0.f);
		
		// get the random sample
		v2f sample = sampler.next2D();

	    // Compute the wi vector based on the sampled vector
		v3f wi = Warp::squareToCosineHemisphere(sample);
		i.wi = wi;

		// Get the pdf value
		*_pdf = pdf(i);

		v3f bsdf = eval(i);

		// compute the Light value
		val = bsdf / *_pdf;

        // TODO(A3): Implement this
        return val;
    }

    std::string toString() const override { return "Diffuse"; }
};

TR_NAMESPACE_END