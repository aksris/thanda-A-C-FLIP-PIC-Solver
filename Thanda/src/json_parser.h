#ifndef JSON_PARSER_H
#define JSON_PARSER_H
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>
#include <json/json.h>
#include "scene/scene.h"
#include <glm/glm.hpp>

class json_parser
{
public:
    json_parser();
    void readFromFile(const char* filename, Scene& s);
};

#endif // JSON_PARSER_H
