#include "json_parser.h"

json_parser::json_parser()
{

}

void json_parser::readFromFile(const char* filename, Scene& s){
    std::string JSON_CONTENT;
    std::ifstream JSONStream(filename, std::ios::in);

    if(JSONStream.is_open()){
        std::string Line = "";
        while(getline(JSONStream, Line))
            JSON_CONTENT += "\n" + Line;
        JSONStream.close();
    }else{
        printf("Impossible to open %s. Are you in the right directory ? Don't forget to read the FAQ !\n", filename);
        getchar();
    }

    Json::Value root;   // will contains the root value after parsing.
    Json::Reader reader;

    bool parsingSuccessful = reader.parse(JSON_CONTENT.c_str(),root);
    if ( !parsingSuccessful )
    {
        // report to the user the failure and their locations in the document.
        std::cout  << "Failed to parse configuration\n"
                   << reader.getFormattedErrorMessages();
        return;
    }

    float scale_x = root["containerDim"]["scaleX"].asFloat();
    float scale_y = root["containerDim"]["scaleY"].asFloat();
    float scale_z = root["containerDim"]["scaleZ"].asFloat();

    float bound_x = root["particleDim"]["boundX"].asFloat();
    float bound_y = root["particleDim"]["boundY"].asFloat();
    float bound_z = root["particleDim"]["boundZ"].asFloat();

    float separation = root["particleSeparation"].asFloat();

    s.containerBounds = glm::vec3(scale_x, scale_y, scale_z);
    s.particleBounds = glm::vec3(bound_x, bound_y, bound_z);
    s.particle_separation = separation;

}

