//
//  viewer.cpp
//  Thanda

#include "viewer.hpp"

Viewer::Viewer(int width, int height, const Scene& s) {
	w = width;
	h = height;
	scene = s;
	fluid = new FluidSolver(scene);
}

Viewer::~Viewer() {
	delete fluid;
}

void Viewer::initializeShader() {

	programID = scene.LoadShaders("src/resources/Particle.vertexshader", "src/resources/Particle.fragmentshader");
	/* for particles and billboards */
	// Vertex shader
	CameraRight_worldspace_ID = glGetUniformLocation(programID, "CameraRight_worldspace");
	CameraUp_worldspace_ID = glGetUniformLocation(programID, "CameraUp_worldspace");
	ViewProjMatrixID = glGetUniformLocation(programID, "VP");
	TextureID = glGetUniformLocation(programID, "myTextureSampler");
	/* for cube */
	programIDGeometry = scene.LoadShaders("src/resources/TransformVertexShader.vertexshader", "src/resources/ColorFragmentShader.fragmentshader");
	// Get a handle for our "MVP" uniform
	/*end cube*/

	/* MVP matrices */
	geomMatrixID = glGetUniformLocation(programIDGeometry, "MVP");

}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode) {
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);

}

void Viewer::initializeGL() {
	if (!glfwInit())
	{
		//fprintf( stderr, "Failed to initialize GLFW\n" );
		getchar();
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef FULLSCRN
	window = glfwCreateWindow(w, h, "Thanda Hai", glfwGetPrimaryMonitor(), NULL);
#else
	window = glfwCreateWindow(w, h, "Thanda Hai", NULL, NULL);
#endif

	if (window == NULL) {
		getchar();
		glfwTerminate();
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		getchar();
		glfwTerminate();
	}

	glfwSetKeyCallback(window, key_callback);
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	// Hide the mouse and enable unlimited mouvement
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	// Set the mouse at the center of the screen
	glfwPollEvents();
	glfwSetCursorPos(window, w / 2, h / 2);

	// Dark blue background
	glClearColor(1.f, 1.f, 1.f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);
}


void Viewer::drawCube() {

	// Enable cubeShader and display cube

	glUseProgram(programIDGeometry);
	glm::mat4 cubeModelMatrix = glm::translate(glm::mat4(1.f), glm::vec3(scene.containerBounds.x / 2.f, scene.containerBounds.y / 2.f, scene.containerBounds.z / 2.f));
	cubeModelMatrix = (glm::scale(cubeModelMatrix, glm::vec3(scene.containerBounds.x-1, scene.containerBounds.y-1, scene.containerBounds.z-1)));
	cube.modelMatrix = cubeModelMatrix;
	glm::mat4 cubeMVP = ProjectionMatrix * ViewMatrix * cubeModelMatrix;

	glUniformMatrix4fv(geomMatrixID, 1, GL_FALSE, &cubeMVP[0][0]);
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, cubevertexbuffer);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);

	glDrawElements(GL_LINES, 36, GL_UNSIGNED_INT, 0);
	glDisableVertexAttribArray(0);


}

void Viewer::drawParticles(int ParticlesCount) {
	/*
	Borrowed from http://www.opengl-tutorial.org/intermediate-tutorials/billboards-particles/particles-instancing/
	*/

	// These functions are specific to glDrawArrays*Instanced*.
	// The first parameter is the attribute buffer we're talking about.
	// The second parameter is the "rate at which generic vertex attributes advance when rendering multiple instances"
	// http://www.opengl.org/sdk/docs/man/xhtml/glVertexAttribDivisor.xml
	glVertexAttribDivisor(0, 0); // particles vertices : always reuse the same 4 vertices -> 0
	glVertexAttribDivisor(1, 1); // positions : one per quad (its center)                 -> 1
	glVertexAttribDivisor(2, 1); // color : one per quad                                  -> 1

	// This is equivalent to :
	// for(i in ParticlesCount) : glDrawArrays(GL_TRIANGLE_STRIP, 0, 4),
	// but faster.
	glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, ParticlesCount);
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(2);
}

void Viewer::allocateParticleBuffers(GLuint billboard_vertex_buffer, GLuint particles_position_buffer, GLuint particles_color_buffer) {
	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
	glVertexAttribPointer(
		0,                  // attribute. No particular reason for 0, but must match the layout in the shader.
		3,                  // size
		GL_FLOAT,           // type
		GL_FALSE,           // normalized?
		0,                  // stride
		(void*)0            // array buffer offset
	);

	// 2nd attribute buffer : positions of particles' centers
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
	glVertexAttribPointer(
		1,                                // attribute. No particular reason for 1, but must match the layout in the shader.
		4,                                // size : x + y + z + size => 4
		GL_FLOAT,                         // type
		GL_FALSE,                         // normalized?
		0,                                // stride
		(void*)0                          // array buffer offset
	);

	// 3rd attribute buffer : particles' colors
	glEnableVertexAttribArray(2);
	glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
	glVertexAttribPointer(
		2,                                // attribute. No particular reason for 1, but must match the layout in the shader.
		4,                                // size : r + g + b + a => 4
		GL_UNSIGNED_BYTE,                 // type
		GL_TRUE,                          // normalized?    *** YES, this means that the unsigned char[4] will be accessible with a vec4 (floats) in the shader ***
		0,                                // stride
		(void*)0                          // array buffer offset
	);
}

void SortParticles(std::vector<Particle>& ParticlesContainer, int MaxParticles)
{
	std::sort(&ParticlesContainer[0], &ParticlesContainer[MaxParticles]);
}

void Viewer::display() {


	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);


	// generate cube and particle buffers once
	glGenBuffers(1, &elementbuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, cube.cub_idx.size() * sizeof(GLuint), &cube.cub_idx[0], GL_STATIC_DRAW);

	glGenBuffers(1, &cubevertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, cubevertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, cube.g_cube_vertex_buffer_data.size() * sizeof(GLfloat), &cube.g_cube_vertex_buffer_data[0], GL_STATIC_DRAW);

	GLuint Texture = scene.loadDDS("src/resources/particle.DDS");
	// The VBO containing the 4 vertices of the particles.
	// Thanks to instancing, they will be shared by all particles.
	static const GLfloat g_vertex_buffer_data[] = {
		 -0.5f, -0.5f, 0.0f,
		  0.5f, -0.5f, 0.0f,
		 -0.5f,  0.5f, 0.0f,
		  0.5f,  0.5f, 0.0f,
	};
	GLuint billboard_vertex_buffer;
	glGenBuffers(1, &billboard_vertex_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, billboard_vertex_buffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data), g_vertex_buffer_data, GL_STATIC_DRAW);

	GLuint particles_position_buffer;
	glGenBuffers(1, &particles_position_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	glBufferData(GL_ARRAY_BUFFER, fluid->MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);

	GLuint particles_color_buffer;
	glGenBuffers(1, &particles_color_buffer);
	glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
	// Initialize with empty (NULL) buffer : it will be updated later, each frame.
	glBufferData(GL_ARRAY_BUFFER, fluid->MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW);


	static GLfloat* g_particule_position_size_data = new GLfloat[fluid->MaxParticles * 4];
	static GLubyte* g_particule_color_data = new GLubyte[fluid->MaxParticles * 4];
	glm::mat4 ViewProjectionMatrix;

	//Step 1 - construct mac grid
	fluid->constructMACGrid();


	// Step 2 - Particle Seeding
	fluid->genParticles();

	string title = "Bounce, baby, bounce: " + std::to_string(fluid->ParticlesContainer.size()) + " Particles";
	glfwSetWindowTitle(window, title.c_str());

	double lastTime = glfwGetTime();
	bool fin = false;
	//update this to have CFL condition
	float dt = scene.step, t = 0.f, frame_time = 1.f / 30.f, step = 0.05f;
	int counter = 0;
	double fps = 0.0;
	bool toggle = false;

	while (!glfwWindowShouldClose(window)) {
		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		//reset after every iteration
		t = 0.f, frame_time = 1.f / 30.f;

		counter++;
		auto currentTime = glfwGetTime();
		auto deltaTime = -(lastTime - currentTime);

		if (deltaTime >= 1.0)
		{
			fps = counter / deltaTime;
			title = "Bounce, baby, bounce: " + std::to_string(fluid->ParticlesContainer.size()) + " Particles | FPS: " + std::to_string(fps) + " | Step: " + std::to_string(dt);
			glfwSetWindowTitle(window, title.c_str());
			counter = 0;
			lastTime = currentTime;
		}

		if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
			toggle = !toggle;
		}

		if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS) {
			// Step 2 - Particle Seeding
			fluid->genParticles();
		}

		if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) {
			// Step 2 - Particle Seeding
			fluid->ParticlesContainer.clear();
			dt = step;
			fluid->genParticles();
		}

		if (glfwGetKey(window, GLFW_KEY_RIGHT_BRACKET) == GLFW_PRESS) {
			// Increase step
			if (dt <= 1.f)
				dt += step;
		}

		if (glfwGetKey(window, GLFW_KEY_LEFT_BRACKET) == GLFW_PRESS) {
			// Decrease step
			if (dt - step > 0.f)
			{
				dt -= 0.05f;
			}
		}

		if (toggle)
		{
			fluid->step(dt);
		}

		// setup camera
		camera.computeMatricesFromInputs(window);
		ViewMatrix = camera.getViewMatrix();
		ProjectionMatrix = camera.getProjectionMatrix();
		CameraPosition = glm::vec3(glm::inverse(ViewMatrix)[3]);

		ViewProjectionMatrix = ProjectionMatrix * ViewMatrix;

		drawCube();

		int ParticlesCount = 0;

		for (int i = 0; i < fluid->ParticlesContainer.size(); i++) {

			Particle& p = fluid->ParticlesContainer[i]; // shortcut
			p.cameradistance = glm::length2(p.pos - CameraPosition);

			// Fill the GPU buffer
			g_particule_position_size_data[4 * ParticlesCount + 0] = p.pos.x;
			g_particule_position_size_data[4 * ParticlesCount + 1] = p.pos.y;
			g_particule_position_size_data[4 * ParticlesCount + 2] = p.pos.z;

			g_particule_position_size_data[4 * ParticlesCount + 3] = p.size;

		vec3 color = glm::normalize(p.speed);
			float maxV = glm::max(glm::max(color.r, color.g), color.b);
			int maxElem = color.r > color.g ? (color.r > color.b ? 0 : 2) : (color.g > color.b ? 1 : 2);
			color = glm::vec3(maxElem != 0 ? 0.f : color.r, maxElem != 1 ? 0.f : color.g, maxElem != 2 ? 0.f : color.b) * 255.f;

			g_particule_color_data[4 * ParticlesCount + 0] = (color.r);
			g_particule_color_data[4 * ParticlesCount + 1] = (color.g);
			g_particule_color_data[4 * ParticlesCount + 2] = (color.b);
			g_particule_color_data[4 * ParticlesCount + 3] = p.a;

			ParticlesCount++;

		}

		int MaxParticles = fluid->ParticlesContainer.size();
		SortParticles(fluid->ParticlesContainer, MaxParticles);

		//  Use our shader
		glUseProgram(programID);

		glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
		glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
		glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLfloat) * 4, g_particule_position_size_data);

		glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
		glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
		glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLubyte) * 4, g_particule_color_data);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		// Bind our texture in Texture Unit 0
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, Texture);
		// Set our "myTextureSampler" sampler to user Texture Unit 0
		glUniform1i(TextureID, 0);

		glUniform3f(CameraRight_worldspace_ID, ViewMatrix[0][0], ViewMatrix[1][0], ViewMatrix[2][0]);
		glUniform3f(CameraUp_worldspace_ID, ViewMatrix[0][1], ViewMatrix[1][1], ViewMatrix[2][1]);
		glUniformMatrix4fv(ViewProjMatrixID, 1, GL_FALSE, &ViewProjectionMatrix[0][0]);

		allocateParticleBuffers(billboard_vertex_buffer, particles_position_buffer, particles_color_buffer);
		drawParticles(ParticlesCount);

		glDisable(GL_BLEND);

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();
	}
	glfwTerminate();
	delete[] g_particule_position_size_data;
	delete[] g_particule_color_data;
}
