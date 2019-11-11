/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Modified Phong reflectance model
 */
struct PhongBSDF : BSDF {

    std::unique_ptr<Texture < v3f>> specularReflectance;
    std::unique_ptr<Texture < v3f>> diffuseReflectance;
    std::unique_ptr<Texture < float>> exponent;
    float specularSamplingWeight;
    float scale;

    PhongBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        const tinyobj::material_t& mat = scene.materials[matID];

        if (mat.specular_texname.empty())
            specularReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.specular)));
        else
            specularReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.specular_texname));

        if (mat.diffuse_texname.empty())
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
        else
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

        exponent = std::unique_ptr<Texture<float>>(new ConstantTexture1f(mat.shininess));

        //get scale value to ensure energy conservation
        v3f maxValue = specularReflectance->getMax() + diffuseReflectance->getMax();
        float actualMax = max(max(maxValue.x, maxValue.y), maxValue.z);
        scale = actualMax > 1.0f ? 0.99f * (1.0f / actualMax) : 1.0f;

        float dAvg = getLuminance(diffuseReflectance->getAverage() * scale);
        float sAvg = getLuminance(specularReflectance->getAverage() * scale);
        specularSamplingWeight = sAvg / (dAvg + sAvg);

        components.push_back(EGlossyReflection);
        components.push_back(EDiffuseReflection);

        combinedType = 0;
        for (unsigned int component : components)
            combinedType |= component;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    v3f eval(const SurfaceInteraction& i) const override {
		v3f val(0.f);
		v3f black_colour(0.f);

		// TODO(A2): Implement this

		// need to normalize both vectors 
		v3f normal_wo = normalize(i.wo);
		v3f normal_wi = normalize(i.wi);

		// Perform to check to see if rays are front facing to the surface.
		if (!(Frame::cosTheta(normal_wo) > 0 && Frame::cosTheta(normal_wi) > 0)) {
			return black_colour;
		}
		
		// get all the values for diffuseReflectance, specularReflectance and exponent 
		v3f rho_d = diffuseReflectance->eval(worldData, i);
		v3f rho_s = specularReflectance->eval(worldData, i);
		float n = exponent->eval(worldData, i);

		// calculate remaining unknowns in Phong BRDF equation
		v3f w_r = PhongBSDF::reflect(normal_wi);// specular reflection
		float cos_alpha = glm::dot(normal_wo, w_r);

		// calculate the phong BRDF function
		v3f phong_BRDF = (rho_d / M_PI) + rho_s * (n + 2) * max(0.f, pow(cos_alpha, n)) / (2 * M_PI) ;
		// multiply phong BRDF by the cos forshotening factor
		val = scale * phong_BRDF * Frame::cosTheta(normal_wi);

        return val;
    }

    float pdf(const SurfaceInteraction& i) const override {
        float pdf = 0.f;
		float n = exponent->eval(worldData, i);

		pdf = Warp::squareToPhongLobePdf(i.wi, n);
        // TODO(A3): Implement this
        return pdf;
    }

    v3f sample(SurfaceInteraction& i, Sampler& sampler, float* _pdf) const override {
        v3f val(0.f);
        // TODO(A3): Implement this

		// get the exponent
		float n = exponent->eval(worldData, i);

		// Get the wi from random sample
		v2f sample = sampler.next2D();
		v3f wi;
		float pdf_d, pdf_p;

		if (sampler.next() <= specularSamplingWeight) {
			wi = Warp::squareToPhongLobe(sample, n);
			i.wi = wi;

			pdf_p = this->pdf(i) * specularSamplingWeight;

			// set new frame so that we can get wi vector in phong lobe
			v3f wr = reflect(i.wo);
			Frame frame_R(wr);
			i.wi = normalize(frame_R.toWorld(wi));

			pdf_d = Warp::squareToCosineHemispherePdf(i.wi) * (1.f - specularSamplingWeight);
		}
		else {
			wi = Warp::squareToCosineHemisphere(sample);
			pdf_d = Warp::squareToCosineHemispherePdf(wi) * (1.f - specularSamplingWeight);
			
			// rotate the around the reflected ray and calculate the phong lobe in pdf format 
			v3f wr = reflect(i.wo);
			Frame frame_R(wr);
			i.wi = normalize(frame_R.toLocal(wi));

			pdf_p = this->pdf(i) * specularSamplingWeight;

			i.wi = wi;
		}		

		*_pdf = pdf_d + pdf_p;

		// ensure pdf is greater than 0
		v3f brdf = eval(i);

		if (*_pdf > 0.f) {
			val = brdf / *_pdf;
		}	

        return val;
    }

    std::string toString() const override { return "Phong"; }
};

TR_NAMESPACE_END