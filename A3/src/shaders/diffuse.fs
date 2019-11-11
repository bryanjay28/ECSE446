/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446 Realistic/Advanced Image Synthesis.
*/

#version 330 core
#define PI       3.14159265358979323846   // pi

// get the desired uniforms
uniform vec3 camPos;
uniform vec3 albedo;
uniform vec3 lightPos;
uniform vec3 lightIntensity;

in vec3 vNormal;
in vec3 vPos;

out vec3 color;

void main()
{
	// get wi
	vec3 wi = lightPos - vPos;
	float dist = distance(lightPos, vPos);
	
	// need to normalize both vectors 
	vec3 normal_n = normalize(vNormal);
	vec3 normal_wi = normalize(wi);

	// do the dot product to get the angle between normal and incident light
	float cos_theta = dot(normal_wi, normal_n);

	// same calculation as in part diffuse.h
	vec3 BSDF = cos_theta * albedo / PI;

	vec3 Li = lightIntensity * BSDF / pow(dist, 2);

	color = Li;
}
