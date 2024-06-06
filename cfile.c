#include "rlEmbree.h"
#include <stdio.h>
int main(int argc, char* argv[]){
    reDrawCube((Vector3){0,0,0}, (Vector3){1,1,1});
    Ray ray = {.position = {3,3,3}, .direction = {-1,-0.99,-0.98}};
    FinishScene();
    struct Rayhit hit = Raycast(ray);
    printf("%f\n", hit.tfar);
    printf("%f, %f, %f\n", hit.position.x, hit.position.y, hit.position.z);
}