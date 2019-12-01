/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <core/core.h>
#include "core/renderpass.h"
#include "tiny_obj_loader.h"
#include "integrators/path.h"

TR_NAMESPACE_BEGIN

/**
 * Global Illumination baking renderpass.
 */
struct GIPass : RenderPass {
    GLuint shader{0};

    GLuint modelMatUniform{0};
    GLuint viewMatUniform{0};
    GLuint projectionMatUniform{0};

    int m_samplePerVertex;

    std::unique_ptr<PathTracerIntegrator> m_ptIntegrator;

    explicit GIPass(const Scene& scene) : RenderPass(scene) {
        m_ptIntegrator = std::unique_ptr<PathTracerIntegrator>(new PathTracerIntegrator(scene));
        m_ptIntegrator->m_maxDepth = scene.config.integratorSettings.gi.maxDepth;
        m_ptIntegrator->m_rrProb = scene.config.integratorSettings.gi.rrProb;
        m_ptIntegrator->m_rrDepth = scene.config.integratorSettings.gi.rrDepth;
        m_samplePerVertex = scene.config.integratorSettings.gi.samplesByVertex;
    }

    virtual void buildVBO(size_t objectIdx) override {
        GLObject& obj = objects[objectIdx];

        // TODO(A5): Implement this

		// set the object vertices in the scene similar to A4 
		obj.nVerts = scene.getObjectNbVertices(objectIdx);
		obj.vertices.resize(N_ATTR_PER_VERT * obj.nVerts);

		// vertex id in array
		int vertex_id = 0;

		// iterate throught the vertices of the object and set it
		for (size_t vertex = 0; vertex < obj.nVerts; vertex++) {
			// create the sampler in the loop to get the same random number pattern
			Sampler sampler = Sampler(260738764);

			v3f normal = scene.getObjectVertexNormal(objectIdx, vertex);
			v3f pos = scene.getObjectVertexPosition(objectIdx, vertex);

			// create the surface interaction object 
			SurfaceInteraction info = SurfaceInteraction();
			
			info.primID = scene.getPrimitiveID(vertex);
			info.matID = scene.getMaterialID(objectIdx, info.primID);
			
			info.frameNg = Frame(normal);
			info.frameNs = Frame(normal);
			info.shapeID = objectIdx;
			info.wo = v3f(0, 0, 1); // set an arbitrary ray direction
			info.p = pos + Epsilon * normal; // shift the shading point to prevent self intersections

			// create the ray to call the render explicit section
			Ray ray(info.p, normal);
			v3f RGB(0.f);

			// sample the vertex according to the given variable and take the average of random samples
			for (int x = 0; x < m_samplePerVertex; x++) {
				RGB += m_ptIntegrator->renderExplicit(ray, sampler, info);
			}

			RGB /= m_samplePerVertex;

			// set the position and colour for the vertex
			obj.vertices[vertex_id] = pos.x;
			obj.vertices[vertex_id + 1] = pos.y;
			obj.vertices[vertex_id + 2] = pos.z;

			obj.vertices[vertex_id + 3] = RGB.x;
			obj.vertices[vertex_id + 4] = RGB.y;
			obj.vertices[vertex_id + 5] = RGB.z;

			vertex_id += N_ATTR_PER_VERT;
		}

        // VBO
        glGenVertexArrays(1, &obj.vao);
        glBindVertexArray(obj.vao);

        glGenBuffers(1, &obj.vbo);
        glBindBuffer(GL_ARRAY_BUFFER, obj.vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     sizeof(GLfloat) * obj.nVerts * N_ATTR_PER_VERT,
                     (GLvoid*) (&obj.vertices[0]),
                     GL_STATIC_DRAW);
    }

    bool init(const Config& config) override {
        RenderPass::init(config);

        // Create shader
        GLuint vs = compileShader("gi.vs", GL_VERTEX_SHADER);
        GLuint fs = compileShader("gi.fs", GL_FRAGMENT_SHADER);
        shader = compileProgram(vs, fs);
        glDeleteShader(vs);
        glDeleteShader(fs);

        // Create uniforms
        modelMatUniform = GLuint(glGetUniformLocation(shader, "model"));
        viewMatUniform = GLuint(glGetUniformLocation(shader, "view"));
        projectionMatUniform = GLuint(glGetUniformLocation(shader, "projection"));

        // Create vertex buffers
        objects.resize(scene.worldData.shapes.size());
        for (size_t i = 0; i < objects.size(); i++) {
            buildVBO(i);
            buildVAO(i);
        }

        return true;
    }

    void cleanUp() override {
        // Delete vertex buffers
        for (size_t i = 0; i < objects.size(); i++) {
            glDeleteBuffers(1, &objects[i].vbo);
            glDeleteVertexArrays(1, &objects[i].vao);
        }

        RenderPass::cleanUp();
    }

    void render() override {
        glBindFramebuffer(GL_FRAMEBUFFER, postprocess_fboScreen);
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);

        // TODO(A5): Implement this

		// similar to normal.h

		// Define shader to use
		glUseProgram(shader);

		// Update camera
		glm::mat4 model, view, projection;
		camera.Update();
		camera.GetMatricies(projection, view, model);

		// Pass uniforms
		glUniformMatrix4fv(modelMatUniform, 1, GL_FALSE, &(modelMat[0][0]));
		glUniformMatrix4fv(viewMatUniform, 1, GL_FALSE, &(view[0][0]));
		glUniformMatrix4fv(projectionMatUniform, 1, GL_FALSE, &(projection[0][0]));

		// Draw
		for (auto& object : objects) {

			// vbo -> vertex buffer object, bind the vertex array
			glBindVertexArray(object.vbo);

			// draw the triangles based on vertices
			glDrawArrays(GL_TRIANGLES, 0, object.nVerts);

			// bind the vertex array to 0 and not the object anymore
			glBindVertexArray(0);
		}

        RenderPass::render();
    }

};

TR_NAMESPACE_END
