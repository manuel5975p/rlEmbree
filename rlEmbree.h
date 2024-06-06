#ifndef RL_EMBREE_H
#define RL_EMBREE_H
#include <raylib.h>
#ifdef __cplusplus
#include <cstddef>
using std::size_t;
extern "C"{
#else
#include <stddef.h>
#endif

struct Rayhit{
    Vector3 position;
    Vector3 normal;
    unsigned int geomID;
    unsigned int faceID;
    float tnear, tfar;
    float u, v;
};

struct sphere_info{
    Vector3 pos;
    float radius;
};
struct triangle_info{
    Vector3 vertices[3];
};
struct cube_info{
    Vector3 from, to;
};
typedef struct RTCGeometryTy* RTCGeometry;
typedef struct RTCDeviceTy* RTCDevice;
struct RTCRayHit;
extern RTCDevice gDevice;
RTCGeometry raylib_to_embree(Mesh mesh);
unsigned int reDrawMesh(Mesh mesh);
unsigned int reDrawCube(Vector3 from, Vector3 end);
unsigned int reDrawCubeTransformed(Vector3 extents, Matrix trf);
unsigned int reDrawSphere(Vector3 pos, float radius);
unsigned int reDrawSpheres(const Vector3* pos, const float* radii, size_t count);
unsigned int reDrawMeshTransformed(Mesh mesh, const Matrix* trf);
unsigned int reDrawTriangle(Vector3 v1, Vector3 v2, Vector3 v3);
unsigned int* reDrawMeshInstanced(Mesh mesh, const Matrix* trf, unsigned int instanceCount);
struct RTCRayHit RaycastEx(Ray r);
struct Rayhit Raycast(Ray r);
void FinishScene();
#ifdef __cplusplus
}
#endif
#endif