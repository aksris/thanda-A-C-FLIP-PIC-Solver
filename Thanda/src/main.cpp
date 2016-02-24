#include "main.h"

using namespace std;


int main()
{
    if( !glfwInit() )
    {
        fprintf( stderr, "Failed to initialize GLFW\n" );
        getchar();
        return -1;
    }

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_RESIZABLE,GL_FALSE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow( 1024, 768, "3D Viewer", NULL, NULL);
    if( window == NULL ){
        fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
        getchar();
        glfwTerminate();
//        return -1;
    }
    glfwMakeContextCurrent(window);

    // Initialize GLEW
    glewExperimental = true; // Needed for core profile
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        getchar();
        glfwTerminate();
        return -1;
    }

    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    // Hide the mouse and enable unlimited mouvement
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // Set the mouse at the center of the screen
    glfwPollEvents();
    glfwSetCursorPos(window, 1024/2, 768/2);

    // Dark blue background
    glClearColor(0.4f, 0.4f, 0.4f, 0.0f);

    // Enable depth test
    glEnable(GL_DEPTH_TEST);
    // Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS);

    //AntTweakBar
//    TwInit(TW_OPENGL_CORE, NULL);

    GLuint programID = LoadShaders( "shaders/Particle.vertexshader", "shaders/Particle.fragmentshader" );
    /* for particles and billboards */
    // Vertex shader
    GLuint CameraRight_worldspace_ID  = glGetUniformLocation(programID, "CameraRight_worldspace");
    GLuint CameraUp_worldspace_ID  = glGetUniformLocation(programID, "CameraUp_worldspace");
    GLuint ViewProjMatrixID = glGetUniformLocation(programID, "VP");

    char* filename = "scene/scene1.json";
    std::vector<float> inputs;
    Scene s;
    json_parser js;
    js.readFromFile(filename, s);

    // fragment shader
    GLuint TextureID  = glGetUniformLocation(programID, "myTextureSampler");
    /*end particles*/
    /* for cube */
    GLuint programIDCube = LoadShaders( "shaders/TransformVertexShader.vertexshader", "shaders/ColorFragmentShader.fragmentshader" );

    // Get a handle for our "MVP" uniform
    /*end cube*/

    /* MVP matrices */
    GLuint cubeMatrixID = glGetUniformLocation(programIDCube, "MVP");
    static const GLfloat g_cube_vertex_buffer_data[] = {
        -0.5f, -0.5f, -0.5f,
        -0.5f, -0.5f, 0.5f,
        -0.5f, 0.5f, 0.5f,
        -0.5f, 0.5f, -0.5f,
        0.5f, -0.5f, -0.5f,
        0.5f, -0.5f, 0.5f,
        0.5f, 0.5f, 0.5f,
        0.5f, 0.5f, -0.5f
    };

    std::vector<GLuint> cub_idx;

    cub_idx.push_back(0);
    cub_idx.push_back(1);
    cub_idx.push_back(1);
    cub_idx.push_back(2);
    cub_idx.push_back(2);
    cub_idx.push_back(3);
    cub_idx.push_back(3);
    cub_idx.push_back(0);

    cub_idx.push_back(0);
    cub_idx.push_back(4);
    cub_idx.push_back(1);
    cub_idx.push_back(5);
    cub_idx.push_back(2);
    cub_idx.push_back(6);
    cub_idx.push_back(3);
    cub_idx.push_back(7);

    cub_idx.push_back(4);
    cub_idx.push_back(5);
    cub_idx.push_back(5);
    cub_idx.push_back(6);
    cub_idx.push_back(6);
    cub_idx.push_back(7);
    cub_idx.push_back(7);
    cub_idx.push_back(4);

    GLuint elementbuffer;
    glGenBuffers(1, &elementbuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, cub_idx.size() * sizeof(GLuint), &cub_idx[0] , GL_STATIC_DRAW);

    GLuint cubevertexbuffer;
    glGenBuffers(1, &cubevertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, cubevertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_cube_vertex_buffer_data), g_cube_vertex_buffer_data, GL_STATIC_DRAW);
    /*end MVP matrices cube*/

    GLuint VertexArrayID;
    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    /*create particles array */
    GLuint VertexArrayIDCube;
    glGenVertexArrays(1, &VertexArrayIDCube);
    glBindVertexArray(VertexArrayIDCube);
    ParticleSystem ps;
    static GLfloat* g_particule_position_size_data = new GLfloat[ps.MaxParticles * 4];
    static GLubyte* g_particule_color_data         = new GLubyte[ps.MaxParticles * 4];

    GLuint Texture = loadDDS("shaders/particle.DDS");
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
    glBufferData(GL_ARRAY_BUFFER, ps.MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW);

    GLuint particles_color_buffer;
    glGenBuffers(1, &particles_color_buffer);
    glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
    // Initialize with empty (NULL) buffer : it will be updated later, each frame.
    glBufferData(GL_ARRAY_BUFFER, ps.MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW);

    //generate particles
    ps.genParticles(s.particle_separation, s.particleBounds.x, s.particleBounds.y, s.particleBounds.z);

    double lastTime = glfwGetTime();
    do {
        // Clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        double currentTime = glfwGetTime();
        double delta = currentTime - lastTime;
        lastTime = currentTime;

        computeMatricesFromInputs(window);
        glm::mat4 ProjectionMatrixParticles = getProjectionMatrix();
        glm::mat4 ViewMatrixParticles = getViewMatrix();
        // We will need the camera's position in order to sort the particles
        // w.r.t the camera's distance.
        // There should be a getCameraPosition() function in common/controls.cpp,
        // but this works too.
        glm::vec3 CameraPosition(glm::inverse(ViewMatrixParticles)[3]);

        glm::mat4 ViewProjectionMatrixParticles = ProjectionMatrixParticles * ViewMatrixParticles;

        //simulation loop
        int ParticlesCount = 0;
        for(int i=0; i< ps.ParticlesContainer.size(); i++){

            Particle& p = ps.ParticlesContainer[i]; // shortcut

            p.speed += glm::vec3(0.0f,-9.81f, 0.0f) * (float)delta * 0.5f;
            p.pos += p.speed * (float)delta;
            p.cameradistance = glm::length2( p.pos - CameraPosition );

            // Fill the GPU buffer
            g_particule_position_size_data[4*ParticlesCount+0] = p.pos.x;
            g_particule_position_size_data[4*ParticlesCount+1] = p.pos.y;
            g_particule_position_size_data[4*ParticlesCount+2] = p.pos.z;

            g_particule_position_size_data[4*ParticlesCount+3] = p.size;

            g_particule_color_data[4*ParticlesCount+0] = p.r;
            g_particule_color_data[4*ParticlesCount+1] = p.g;
            g_particule_color_data[4*ParticlesCount+2] = p.b;
            g_particule_color_data[4*ParticlesCount+3] = p.a;

            ParticlesCount++;

        }
//        ps.sortParticles();
       int MaxParticles = ps.ParticlesContainer.size();
       glBindBuffer(GL_ARRAY_BUFFER, particles_position_buffer);
       glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLfloat), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
       glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLfloat) * 4, g_particule_position_size_data);

       glBindBuffer(GL_ARRAY_BUFFER, particles_color_buffer);
       glBufferData(GL_ARRAY_BUFFER, MaxParticles * 4 * sizeof(GLubyte), NULL, GL_STREAM_DRAW); // Buffer orphaning, a common way to improve streaming perf. See above link for details.
       glBufferSubData(GL_ARRAY_BUFFER, 0, ParticlesCount * sizeof(GLubyte) * 4, g_particule_color_data);

       glEnable(GL_BLEND);
       glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

       // Use our shader
       glUseProgram(programID);

       // Bind our texture in Texture Unit 0
       glActiveTexture(GL_TEXTURE0);
       glBindTexture(GL_TEXTURE_2D, Texture);
       // Set our "myTextureSampler" sampler to user Texture Unit 0
       glUniform1i(TextureID, 0);

       // Same as the billboards tutorial
       glUniform3f(CameraRight_worldspace_ID, ViewMatrixParticles[0][0], ViewMatrixParticles[1][0], ViewMatrixParticles[2][0]);
       glUniform3f(CameraUp_worldspace_ID   , ViewMatrixParticles[0][1], ViewMatrixParticles[1][1], ViewMatrixParticles[2][1]);

       glUniformMatrix4fv(ViewProjMatrixID, 1, GL_FALSE, &ViewProjectionMatrixParticles[0][0]);

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

        //cube
        glDisable(GL_BLEND);
        glUseProgram(programIDCube);
        glm::mat4 cubeModelMatrix(glm::scale(glm::mat4(1.f),glm::vec3(s.containerBounds.x, s.containerBounds.y, s.containerBounds.z )));
        glm::mat4 cubeMVP = ProjectionMatrixParticles * ViewMatrixParticles * cubeModelMatrix;
        glUniformMatrix4fv(cubeMatrixID, 1, GL_FALSE, &cubeMVP[0][0]);
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, cubevertexbuffer);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0  );

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);

        glDrawElements(GL_LINES, 36, GL_UNSIGNED_INT, 0);
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);


        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();

    }while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
            glfwWindowShouldClose(window) == 0 );

    delete[] g_particule_position_size_data;

    // Cleanup VBO and shader
    glDeleteBuffers(1, &particles_color_buffer);
    glDeleteBuffers(1, &particles_position_buffer);
    glDeleteBuffers(1, &billboard_vertex_buffer);
    glDeleteProgram(programID);
    glDeleteTextures(1, &TextureID);
    glDeleteVertexArrays(1, &VertexArrayID);

    // Cleanup VBO and shader
    glDeleteBuffers(1, &cubevertexbuffer);
    glDeleteBuffers(1, &elementbuffer);
    glDeleteProgram(programIDCube);
    glDeleteVertexArrays(1, &VertexArrayIDCube);

    // Close OpenGL window and terminate GLFW
    glfwTerminate();
    return 0;
}

