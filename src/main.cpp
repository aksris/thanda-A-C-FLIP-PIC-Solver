#include "main.hpp"

using namespace std;


int main()
{
    Scene scene;
    char* filename = "src/resources/scene.json";

    scene.parseScene(filename);
    char* obj = "src/resources/cube.obj";
    //scene.LoadOBJ(obj, scene);
    int width = 1024;
    int height = 768;
    Viewer viewer(width, height, scene);

    viewer.initializeGL();
    viewer.initializeShader();
    viewer.display();

    return 0;
}
