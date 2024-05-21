#include <numeric>
#include <raylib.h>
#include <raymath.h>
#include <rlgl.h>
#include <cstring>
#include <vector>
#include <optional>
#include <unordered_map>
#include <iostream>
#include <limits>
#include <cassert>
#include <embree4/rtcore.h>
RTCDevice gDevice;
struct reScene{
    RTCScene m_scene;
    operator RTCScene()const noexcept{
        return m_scene;
    }
    /// Maps geometry IDs to meshes
    std::unordered_map<unsigned int, std::pair<Mesh, std::optional<Matrix>>> geometry_mesh_map;
    void clear(){
        rtcReleaseScene(m_scene);
        m_scene = rtcNewScene(gDevice);
        geometry_mesh_map.clear();
    }
};
reScene gScene;

RTCGeometry raylib_to_embree(Mesh mesh){
    RTCGeometry ret = rtcNewGeometry(gDevice, RTC_GEOMETRY_TYPE_TRIANGLE);
    float* vb = (float*) rtcSetNewGeometryBuffer(ret, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,  3 * sizeof(float), mesh.vertexCount);

    std::memcpy(vb, mesh.vertices, sizeof(float) * mesh.vertexCount * 3);

    if(mesh.indices != nullptr){
        abort();
        unsigned int* ib = (unsigned int*) rtcSetNewGeometryBuffer(ret, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,  3 * sizeof(unsigned), mesh.triangleCount);
        std::memcpy(ib, mesh.vertices, sizeof(unsigned) * mesh.triangleCount * 3);
    }
    else{
        assert(false && "Indexed meshes not supported");
        unsigned int* ib = (unsigned int*) rtcSetNewGeometryBuffer(ret, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,  3 * sizeof(unsigned), mesh.vertexCount);
        std::iota(ib, ib + mesh.vertexCount * 3, 0u);
    }
    
    //rtcSetGeometryTransform(ret, 1, RTC_FORMAT_FLOAT3X4_ROW_MAJOR, trf);
    rtcCommitGeometry(ret);

    return ret;
}
unsigned int reDrawMesh(Mesh mesh){
    RTCGeometry geom = raylib_to_embree(mesh);
    unsigned int ret = rtcAttachGeometry(gScene, geom);
    gScene.geometry_mesh_map[ret] = {mesh, std::nullopt};
    rtcReleaseGeometry(geom);
    return ret;
}
unsigned int reDrawMeshTransformed(Mesh mesh, const Matrix& trf){
    RTCGeometry geom = raylib_to_embree(mesh);
    RTCScene iScene = rtcNewScene(gDevice);
    rtcAttachGeometry(iScene, geom);
    rtcCommitScene(iScene);
    rtcReleaseGeometry(geom);
    RTCGeometry igeom = rtcNewGeometry(gDevice, RTC_GEOMETRY_TYPE_INSTANCE);
    rtcSetGeometryInstancedScene(igeom, iScene);
    rtcSetGeometryTransform(igeom, 0, RTC_FORMAT_FLOAT3X4_ROW_MAJOR, reinterpret_cast<const float*>(&trf));
    rtcCommitGeometry(igeom);
    unsigned int ret = rtcAttachGeometry(gScene, igeom);
    rtcReleaseGeometry(igeom);
    rtcReleaseScene(iScene);
    gScene.geometry_mesh_map[ret] = {mesh, trf};
    return ret;
}

std::vector<unsigned int> reDrawMeshInstanced(Mesh mesh, const Matrix* trf, unsigned int instanceCount){
    RTCGeometry geom = raylib_to_embree(mesh);
    RTCScene iScene = rtcNewScene(gDevice);
    rtcAttachGeometry(iScene, geom);
    rtcCommitScene(iScene);
    rtcReleaseGeometry(geom);
    std::vector<unsigned int> ret(instanceCount);
    for(unsigned int i = 0;i < instanceCount;i++){
        RTCGeometry igeom = rtcNewGeometry(gDevice, RTC_GEOMETRY_TYPE_INSTANCE);
        rtcSetGeometryInstancedScene(igeom, iScene);
        rtcSetGeometryTransform(igeom, 0, RTC_FORMAT_FLOAT3X4_ROW_MAJOR, reinterpret_cast<const float*>(trf + i));
        rtcCommitGeometry(igeom);
        ret[i] = rtcAttachGeometry(gScene, igeom);
        gScene.geometry_mesh_map[ret[i]] = {mesh, trf[i]};
        rtcReleaseGeometry(igeom);
    }
    rtcReleaseScene(iScene);
    return ret;
}
void ClearScene(){
    gScene.clear();
}
struct Rayhit{
    Vector3 position;
    Vector3 normal;
    float tnear, tfar;
    float u, v;
};
RTCRayHit RaycastEx(Ray r){
    RTCRayHit e_rayhit; 
    e_rayhit.ray.org_x = r.position.x; e_rayhit.ray.org_y = r.position.y; e_rayhit.ray.org_z = r.position.z;
    e_rayhit.ray.dir_x = r.direction.x; e_rayhit.ray.dir_y = r.direction.y; e_rayhit.ray.dir_z = r.direction.z;
    
    e_rayhit.ray.tnear = 0.0f;
    e_rayhit.ray.tfar = std::numeric_limits<float>::infinity();
    e_rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
    e_rayhit.ray.mask = -1;
    e_rayhit.ray.flags = 0;
    e_rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    //RTCIntersectArguments context;
    //rtcInitIntersectArguments(&context);
    rtcIntersect1(gScene, &e_rayhit);
    return e_rayhit;
}
inline Vector3 Vector3Transform_nt(Vector3 v, Matrix mat)
{
    Vector3 result = { 0 };

    float x = v.x;
    float y = v.y;
    float z = v.z;

    result.x = mat.m0*x + mat.m4*y + mat.m8*z ;
    result.y = mat.m1*x + mat.m5*y + mat.m9*z ;
    result.z = mat.m2*x + mat.m6*y + mat.m10*z;

    return result;
}
Rayhit Raycast(Ray r){
    RTCRayHit hitex = RaycastEx(r);
    auto it = gScene.geometry_mesh_map.find(hitex.hit.instID[0]);
    Rayhit ret{};
    ret.tnear = hitex.ray.tnear;
    ret.tfar = hitex.ray.tfar;
    if(it != gScene.geometry_mesh_map.end()){
        unsigned int faceID = hitex.hit.primID;
        Vector2 uvv1{it->second.first.texcoords[(faceID * 3 + 0) * 2], it->second.first.texcoords[(faceID * 3 + 0) * 2 + 1]};
        Vector2 uvv2{it->second.first.texcoords[(faceID * 3 + 1) * 2], it->second.first.texcoords[(faceID * 3 + 1) * 2 + 1]};
        Vector2 uvv3{it->second.first.texcoords[(faceID * 3 + 2) * 2], it->second.first.texcoords[(faceID * 3 + 2) * 2 + 1]};
        Vector3 nrv1{it->second.first.normals  [(faceID * 3 + 0) * 3], it->second.first.normals  [(faceID * 3 + 0) * 3 + 1], it->second.first.normals [(faceID * 3 + 0) * 3 + 2]};
        Vector3 nrv2{it->second.first.normals  [(faceID * 3 + 1) * 3], it->second.first.normals  [(faceID * 3 + 1) * 3 + 1], it->second.first.normals [(faceID * 3 + 1) * 3 + 2]};
        Vector3 nrv3{it->second.first.normals  [(faceID * 3 + 2) * 3], it->second.first.normals  [(faceID * 3 + 2) * 3 + 1], it->second.first.normals [(faceID * 3 + 2) * 3 + 2]};
        Vector3 pv1{it->second.first.vertices  [(faceID * 3 + 0) * 3], it->second.first.vertices [(faceID * 3 + 0) * 3 + 1], it->second.first.vertices[(faceID * 3 + 0) * 3 + 2]};
        Vector3 pv2{it->second.first.vertices  [(faceID * 3 + 1) * 3], it->second.first.vertices [(faceID * 3 + 1) * 3 + 1], it->second.first.vertices[(faceID * 3 + 1) * 3 + 2]};
        Vector3 pv3{it->second.first.vertices  [(faceID * 3 + 2) * 3], it->second.first.vertices [(faceID * 3 + 2) * 3 + 1], it->second.first.vertices[(faceID * 3 + 2) * 3 + 2]};
        float x = hitex.hit.u;
        float y = hitex.hit.v;
        float z = 1.0f - x - y;
        Vector2 uv = Vector2Add(Vector2Scale(uvv1, x), Vector2Add(Vector2Scale(uvv2, y), Vector2Scale(uvv3, z)));
        Vector3 nrm = Vector3Add(Vector3Scale(nrv1, x), Vector3Add(Vector3Scale(nrv2, y), Vector3Scale(nrv3, z)));
        Vector3 pm = Vector3Add(Vector3Scale(pv1, x), Vector3Add(Vector3Scale(pv2, y), Vector3Scale(pv3, z)));
        if(it->second.second){
            nrm = Vector3Transform_nt(nrm, it->second.second.value());
            pm  = Vector3Transform   (pm , it->second.second.value());
        }
        ret.u = uv.x;
        ret.v = uv.y;
        ret.position = pm;//Vector3Add(r.position, Vector3Scale(r.direction, hitex.ray.tfar));
        ret.normal = nrm;
    }
    else{
        //abort();
    }
    return ret;
}

int main(){
    
    gDevice = rtcNewDevice(NULL);
    gScene.m_scene = rtcNewScene(gDevice);
    
    InitWindow(1920, 1080, "");
    SetTargetFPS(60);
    Model cube = LoadModel("cube.obj");
    Model human = LoadModel("human.obj");
    Model cottage = LoadModel("cottage.obj");
    rlDisableBackfaceCulling();
    //RTCGeometry geom = raylib_to_embree(mod.meshes[0]);
    Matrix trfs[4] = {
        MatrixTranslate(0, 0, 0.0),
        MatrixTranslate(0, 0, 0.2),
        MatrixTranslate(0, 0, 0.2),
        MatrixTranslate(0, 0,   2),
    };
    //reDrawMeshTransformed(human.meshes[0], MatrixTranslate(1.3765, -15.5177, 1.7395));
    //for(int i = 0;i < cottage.meshCount;i++)
    //    reDrawMeshTransformed(cottage.meshes[i], MatrixTranslate(18.569975, -0.665788, 0.528871));
    //reDrawMeshTransformed(mod.meshes[0], MatrixTranslate(10, 10, 2));
    //reDrawMeshTransformed(mod.meshes[0], MatrixTranslate(10, 10, 2));
    //reDrawMeshTransformed(mod.meshes[0], MatrixTranslate(10, 0, 2));
    //auto vec = reDrawMeshInstanced(mod.meshes[0], trfs, 4);
    
    reDrawMeshTransformed(cube.meshes[0], MatrixTranslate(0,0,0));
    //reDrawMeshTransformed(mod.meshes[0], trfs[0]);
    rtcCommitScene(gScene);
    Camera3D cam;
    cam.fovy = 45.0f;
    cam.projection = CAMERA_PERSPECTIVE;
    cam.position = Vector3{-4,2,-4};
    cam.target = Vector3{0,0,0};
    cam.up = Vector3{0,1,0};
    Image rtImage = GenImageColor(GetScreenWidth(), GetScreenHeight(), BLACK);
    Texture tex = LoadTextureFromImage(rtImage);
    
    while(!WindowShouldClose()){
        ClearBackground(Color { 50, 50, 70, 255 });
        BeginMode3D(cam);
        DrawMesh(cube.meshes[0], LoadMaterialDefault(), MatrixTranslate(0,0,0));
        EndMode3D();
        DrawText(TextFormat("FPS: %d", GetFPS()), 0, 0, 50, GREEN);
        Ray r{.position  = Vector3{0, 0.0, 1000},
              .direction = Vector3{0.0, 0.0, -1}
            };
        RTCRayHit hitex = RaycastEx(r);
        DrawText(TextFormat("hitex.ray.tfar  : %f", hitex.ray.tfar  ), 0, 100, 40, YELLOW);
        DrawText(TextFormat("hitex.hit.geomID: %d", hitex.hit.geomID), 0, 150, 40, YELLOW);
        DrawText(TextFormat("hitex.hit.instID[0]: %d", hitex.hit.instID[0]), 0, 200, 40, YELLOW);
        DrawText(TextFormat("hitex.hit.primID: %d", hitex.hit.primID), 0, 250, 40, YELLOW);
        DrawText(TextFormat("hitex.hit.instPrimID[0]: %d", hitex.hit.instPrimID[0]), 0, 300, 40, YELLOW);
        Rayhit hit = Raycast(r);
        DrawText(TextFormat("hit.tfar: %f", hit.tfar), 600, 100, 40, YELLOW);
        DrawText(TextFormat("hit.uv   : %f, %f", hit.u, hit.v), 600, 150, 40, YELLOW);
        DrawText(TextFormat("hit.pos   : (%f, %f, %f)", hit.position.x, hit.position.y, hit.position.z), 600, 200, 40, YELLOW);
        DrawText(TextFormat("hit.posacc: (%f, %f, %f)", hit.position.x, hit.position.y, hit.position.z), 600, 250, 40, YELLOW);
        DrawText(TextFormat("hit.normal: (%f, %f, %f)", hit.normal.x, hit.normal.y, hit.normal.z), 600, 300, 40, YELLOW);
        //DrawText(TextFormat("hit.hit.instID[0]: %d", hitex.hit.instID[0]),        600, 200, 40, YELLOW);
        //DrawText(TextFormat("hit.hit.primID: %d", hitex.hit.primID),              600, 250, 40, YELLOW);
        //DrawText(TextFormat("hit.hit.instPrimID[0]: %d", hitex.hit.instPrimID[0]),600, 300, 40, YELLOW);
        UpdateTexture(tex, rtImage.data);
        DrawTexture(tex,0,0,WHITE);
        EndDrawing();
    }
}