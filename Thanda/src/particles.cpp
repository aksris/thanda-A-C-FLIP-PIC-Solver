#include "particles.h"

Particle::Particle(){
    pos = glm::vec3(0.f, 0.f, 0.f);
    speed = glm::vec3(0.f, 0.f, 0.f);
    r = 0;
    g = 0;
    b = 220;
    a = 230;
    size = 0.1f;
    angle = 45.f;
    weight = 0.f;
    life = 1.f;
    cameradistance = 10.f;
	density = 0;
	pressure = 0;
}

ParticleSystem::ParticleSystem(){
    LastUsedParticle = 0;
    MaxParticles = 200000;
}

int ParticleSystem::findUnusedParticles(){

    for(int i=LastUsedParticle; i<MaxParticles; i++){
        if (ParticlesContainer[i].life < 0){
            LastUsedParticle = i;
            return i;
        }
    }

    for(int i=0; i<LastUsedParticle; i++){
        if (ParticlesContainer[i].life < 0){
            LastUsedParticle = i;
            return i;
        }
    }

    return 0; // All particles are taken, override the first one
}

void ParticleSystem::sortParticles(){
    std::sort(&ParticlesContainer[0], &ParticlesContainer[MaxParticles]);
}

void ParticleSystem::particlesInit(){
    for(int i=0; i<MaxParticles; i++){
//        ParticlesContainer[i].life = -1.0f;
        ParticlesContainer[i].cameradistance = -1.0f;
        ParticlesContainer[i].size = 0.1f;
    }
}

void ParticleSystem::genParticles(float particle_separation, float boundx, float boundy, float boundz){
    Particle p;
    for(float i = 0; i < boundx; i+= particle_separation){
        for(float j = 0; j < boundy ; j+= particle_separation){
            for(float k = 0; k <boundz ; k+= particle_separation){
                p.pos = glm::vec3(i, j, k);
                ParticlesContainer.push_back(p);
            }
        }
    }
}
//void ParticleSystem::simulateParticles(std::vector<GLFloat>& particle_size_data,std::vector<GLubyte>& particle_color_data, glm::vec3 CameraPosition, double delta){
//    int ParticlesCount = 0;
//    particle_size_data.resize(MaxParticles * 4);
//    particle_color_data.resize(MaxParticles * 4);

//    for(int i=0; i<MaxParticles; i++){

//        Particle& p = ParticlesContainer[i]; // shortcut

//        if(p.life > 0.0f){

//            // Decrease life
//            p.life -= delta;
//            if (p.life > 0.0f){

//                // Simulate simple physics : gravity only, no collisions
//                p.speed += glm::vec3(0.0f,-9.81f, 0.0f) * (float)delta * 0.5f;
//                p.pos += p.speed * (float)delta;
//                p.cameradistance = glm::length2( p.pos - CameraPosition );
//                //ParticlesContainer[i].pos += glm::vec3(0.0f,10.0f, 0.0f) * (float)delta;

//                // Fill the GPU buffer
//                g_particule_position_size_data[4*ParticlesCount+0] = p.pos.x;
//                g_particule_position_size_data[4*ParticlesCount+1] = p.pos.y;
//                g_particule_position_size_data[4*ParticlesCount+2] = p.pos.z;

//                g_particule_position_size_data[4*ParticlesCount+3] = p.size;

//                g_particule_color_data[4*ParticlesCount+0] = p.r;
//                g_particule_color_data[4*ParticlesCount+1] = p.g;
//                g_particule_color_data[4*ParticlesCount+2] = p.b;
//                g_particule_color_data[4*ParticlesCount+3] = p.a;

//            }else{
//                // Particles that just died will be put at the end of the buffer in SortParticles();
//                p.cameradistance = -1.0f;
//            }

//            ParticlesCount++;

//        }
//    }

//    this->sortParticles();
//}
