/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446 Realistic/Advanced Image Synthesis.
*/

#version 330 core
#define PI       3.14159265358979323846   // pi

// get the desired uniforms
uniform vec3 camPos;
uniform vec3 rho_d;
uniform vec3 rho_s;
uniform float exponent;
uniform vec3 lightPos;
uniform vec3 lightIntensity;

in vec3 vNormal;
in vec3 vPos;

out vec3 color;

void main()
{
	// get wi
	vec3 wo = camPos - vPos;
	vec3 wi = lightPos - vPos;
	float dist = distance(lightPos, vPos);
	
	// need to normalize both vectors 
	vec3 normal_n = normalize(vNormal);
	vec3 normal_wo = normalize(wo);
	vec3 normal_wi = normalize(wi);

	// do the dot product to get the angle between normal and incident light
	float cos_theta = max(0.0, dot(normal_n, normal_wi));

	// get cos_alpha with dot product and the relected wave
	vec3 w_r = reflect(-normal_wi, normal_n);
	float cos_alpha = dot(normalize(w_r), normal_wo);


	vec3 phong_BRDF = cos_theta * ((rho_d / PI) + rho_s * (exponent + 2) * max(0.f, pow(cos_alpha, exponent)) / (2 * PI));

	vec3 Li = lightIntensity * phong_BRDF / pow(dist, 2);

	color = Li;
}
