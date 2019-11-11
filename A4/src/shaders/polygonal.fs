/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#version 330 core


#define PI 3.14159265359
#define MAX_NUM_EMITTER_TRIANGLES 40 // Max number of emitter triangles allowed (tuned for A4)
uniform float emitterVertices[MAX_NUM_EMITTER_TRIANGLES*3*3]; // Need to specify a size at compile time (max size is 512 floats)

uniform int nbTriangles; // Use nbTriangles to index the buffer correctly
uniform vec3 lightIrradiance;
uniform vec3 albedo;
uniform vec2 windowSize; // [width, height]

uniform sampler2D cvTerm; // Creates a 2D texture we can sample from (range x,y = [0,1])

in vec3 vNormal;
in vec3 vPos;

out vec3 color;

// Compute edge (v1--v2) contribution
float getEdgeContrib(vec3 v1, vec3 v2, vec3 pos) {
	// Adapt your getEdgeContrib code from the offline part
	float value = 0.f;
	// TODO(A4): Implement this

	vec3 temp1 = normalize(v1 - pos);
	vec3 temp2 = normalize(v2 - pos);
		
	vec3 gamma_i =  normalize(cross(temp1, temp2));
	float theta = acos(dot(temp1, temp2));
		
	value = theta * dot(-gamma_i, vNormal);

	return value;
}


void main()
{	
	// 1) Extract vertices of triangles from `emitterVertices` buffer using `nbTriangles`
	// 2) Calculate G term
	// 3) Subtract modification term for G after extracting it from texture (use built-in `texture()` function)
	//	    e.g. `vec3 delta = texture(cvTerm, coords).xyz;`

	color = vec3(0);
    // TODO(A4): Implement this
	float sum = 0.f;

	for (int i = 0; i < nbTriangles; i++) {
		int id = i * 9;
		vec3 t0 = vec3(emitterVertices[id], emitterVertices[id + 1], emitterVertices[id + 2]);
		vec3 t1 = vec3(emitterVertices[id + 3], emitterVertices[id + 4], emitterVertices[id + 5]);
		vec3 t2 = vec3(emitterVertices[id + 6], emitterVertices[id + 7], emitterVertices[id + 8]);
		sum += getEdgeContrib(t0, t1, vPos);
		sum += getEdgeContrib(t1, t2, vPos);
		sum += getEdgeContrib(t2, t0, vPos);
	}

	vec3 E_poly = lightIrradiance * sum / (2.f * PI);
	vec3 G = albedo * E_poly / PI;

	vec2 coord = vec2(0);

	coord.x = gl_FragCoord.x / windowSize.x;
	coord.y = gl_FragCoord.y / windowSize.y;

	vec3 delta = texture(cvTerm, coord).xyz;

	color = max(vec3(0, 0, 0), G - delta);
}

