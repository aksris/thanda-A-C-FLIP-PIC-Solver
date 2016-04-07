#include "main.hpp"

using namespace std;


int main()
{
    Scene scene;
    char* filename = "src/resources/scene.json";

    scene.parseScene(filename, scene);
    int width = 1920;
    int height = 1080;
    Viewer viewer(width, height, scene);

    viewer.initializeGL();
    viewer.initializeShader();
    viewer.display();

    return 0;
}

