#include <numeric>
#include <raylib.h>
#include <raymath.h>
#include <cstdint>
#include <rlgl.h>
#include <cstring>
#include <vector>
#include <algorithm>
#include <format>
#include <variant>
#include <optional>
#include <unordered_map>
#include <iostream>
#include <limits>
#include <cassert>
#include <embree4/rtcore.h>
#include "rm_operator.hpp"
RTCDevice gDevice;

struct sphere_info{
    Vector3 pos;
    float radius;
};
struct triangle_info{

};
template <class... Fs> struct Overload : Fs... { using Fs::operator()...; };
template <class... Fs> Overload(Fs...) -> Overload<Fs...>;

struct sceneEntry{
    std::variant<Mesh, sphere_info, triangle_info> geometry;
    std::optional<Matrix> transform;
};
struct reScene{
    RTCScene m_scene;
    operator RTCScene()const noexcept{
        return m_scene;
    }
    /// Maps geometry IDs to meshes
    std::unordered_map<unsigned int, sceneEntry> geometry_mesh_map;
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
unsigned int reDrawSphere(Vector3 pos, float radius){
    RTCGeometry geom = rtcNewGeometry(gDevice, RTC_GEOMETRY_TYPE_SPHERE_POINT);
    float* vb = (float*) rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4,  4 * sizeof(float), 1);
    vb[0] = pos.x;
    vb[1] = pos.y;
    vb[2] = pos.z;
    vb[3] = radius;
    rtcCommitGeometry(geom);
    unsigned int ret = rtcAttachGeometry(gScene, geom);
    
    gScene.geometry_mesh_map[ret] = sceneEntry{sphere_info{pos, radius}, std::nullopt};
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
unsigned int reDrawTriangle(Vector3 v1, Vector3 v2, Vector3 v3){
    RTCGeometry geom = rtcNewGeometry(gDevice, RTC_GEOMETRY_TYPE_TRIANGLE);
    float* vb = (float*) rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,  3 * sizeof(float), 3);
    //float* vb2 = (float*) rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_NORMAL, 0, RTC_FORMAT_FLOAT3,  3 * sizeof(float), 3);
    Vector3 vs[3] = {v1, v2, v3};
    //Vector3 normal = Vector3CrossProduct(v2 - v1, v3 - v1);
    std::memcpy(vb, vs, sizeof(float) * 3 * 3);
    unsigned int* ib = (unsigned int*) rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,  3 * sizeof(unsigned), 3);
    std::iota(ib, ib + 3 * 3, 0u);
    rtcCommitGeometry(geom);
    auto geoid = rtcAttachGeometry(gScene, geom);
    gScene.geometry_mesh_map[geoid] = sceneEntry{.geometry = triangle_info{}, .transform = std::nullopt};
    rtcReleaseGeometry(geom);
    return geoid;

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
        //std::cout << ret[i] << "\n";
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
    auto it = gScene.geometry_mesh_map.find(hitex.hit.geomID + (hitex.hit.instID[0] == RTC_INVALID_GEOMETRY_ID ? 0 : hitex.hit.instID[0]));
    Rayhit ret{};
    ret.tnear = hitex.ray.tnear;
    ret.tfar = hitex.ray.tfar;
    //if(!std::isinf(hitex.ray.tfar)){
    //    std::cout << "Geomid: " << hitex.hit.geomID << "\n";
    //    std::cout << "Instanceid: " << hitex.hit.instID[0] << "\n";
    //}
    //std::cout << r.direction << "\n";
    if(it != gScene.geometry_mesh_map.end()){
        
        unsigned int faceID = hitex.hit.primID;
        std::visit(Overload{
            [&](triangle_info sinfo){
                Vector3 p = r.position + r.direction * ret.tfar;
                Vector3 normal{hitex.hit.Ng_x, hitex.hit.Ng_y, hitex.hit.Ng_z};
                ret.position = p;
                ret.normal = Vector3Normalize(normal);
                //abort();
            },
            [&](sphere_info sinfo){
                Vector3 p = r.position + r.direction * ret.tfar;
                Vector3 normal = Vector3Normalize(p - sinfo.pos);
                ret.position = p;
                ret.normal = normal;
                //abort();
            },
            [&](Mesh mesh){
                Vector2 uvv1{mesh.texcoords[(faceID * 3 + 0) * 2], mesh.texcoords[(faceID * 3 + 0) * 2 + 1]};
                Vector2 uvv2{mesh.texcoords[(faceID * 3 + 1) * 2], mesh.texcoords[(faceID * 3 + 1) * 2 + 1]};
                Vector2 uvv3{mesh.texcoords[(faceID * 3 + 2) * 2], mesh.texcoords[(faceID * 3 + 2) * 2 + 1]};
                Vector3 nrv1{mesh.normals  [(faceID * 3 + 0) * 3], mesh.normals  [(faceID * 3 + 0) * 3 + 1], mesh.normals [(faceID * 3 + 0) * 3 + 2]};
                Vector3 nrv2{mesh.normals  [(faceID * 3 + 1) * 3], mesh.normals  [(faceID * 3 + 1) * 3 + 1], mesh.normals [(faceID * 3 + 1) * 3 + 2]};
                Vector3 nrv3{mesh.normals  [(faceID * 3 + 2) * 3], mesh.normals  [(faceID * 3 + 2) * 3 + 1], mesh.normals [(faceID * 3 + 2) * 3 + 2]};
                Vector3 pv1 {mesh.vertices [(faceID * 3 + 0) * 3], mesh.vertices [(faceID * 3 + 0) * 3 + 1], mesh.vertices[(faceID * 3 + 0) * 3 + 2]};
                Vector3 pv2 {mesh.vertices [(faceID * 3 + 1) * 3], mesh.vertices [(faceID * 3 + 1) * 3 + 1], mesh.vertices[(faceID * 3 + 1) * 3 + 2]};
                Vector3 pv3 {mesh.vertices [(faceID * 3 + 2) * 3], mesh.vertices [(faceID * 3 + 2) * 3 + 1], mesh.vertices[(faceID * 3 + 2) * 3 + 2]};
                float x = hitex.hit.u;
                float y = hitex.hit.v;
                float z = 1.0f - x - y;
                //std::cout << std::format("{} {}\n", x, y);
                Vector2 uv = Vector2Add(Vector2Scale(uvv1, x), Vector2Add(Vector2Scale(uvv2, y), Vector2Scale(uvv3, z)));
                Vector3 nrm = Vector3Add(Vector3Scale(nrv1, x), Vector3Add(Vector3Scale(nrv2, y), Vector3Scale(nrv3, z)));
                Vector3 pm = Vector3Add(Vector3Scale(pv1, x), Vector3Add(Vector3Scale(pv2, y), Vector3Scale(pv3, z)));
                if(it->second.transform){
                    nrm = Vector3Transform_nt(nrm, it->second.transform.value());
                    pm  = Vector3Transform   (pm , it->second.transform.value());
                }
                ret.u = uv.x;
                ret.v = uv.y;
                //std::cout << std::format("{}, {}, {}\n", x, y, z);
                //ret.normal = nrm;
                ret.normal = Vector3Transform_nt(Vector3{hitex.hit.Ng_x, hitex.hit.Ng_y, hitex.hit.Ng_z}, it->second.transform.value());
                ret.position = pm;//Vector3Add(r.position, Vector3Scale(r.direction, hitex.ray.tfar));
            },
        },it->second.geometry);
        
    }
    else{
        //abort();
    }
    return ret;
}

void screenloop(const auto& c, int width, int height){
    #pragma omp parallel for collapse(2)
    for(int i = 0;i < height;i++){
        for(int j = 0;j < width;j++){
            c(i,j);
        }
    }
}
int main(){
    
    gDevice = rtcNewDevice(NULL);
    gScene.m_scene = rtcNewScene(gDevice);
    rtcSetSceneBuildQuality(gScene.m_scene, RTC_BUILD_QUALITY_HIGH);
    InitWindow(2560, 1440, "Raylib Raytracer");
    Model cube = LoadModel("cube.obj");
    Model human = LoadModel("ajax.obj");
    Model cottage = LoadModel("cottage.obj");
    rlDisableBackfaceCulling();
    //RTCGeometry geom = raylib_to_embree(mod.meshes[0]);
    Matrix trfs[4] = {
        MatrixTranslate(0, 0, 0.0),
        MatrixTranslate(0, 0, 0.2),
        MatrixTranslate(0, 0, 0.2),
        MatrixTranslate(0, 0,   2),
    };
    size_t meshcount = 10;
    Matrix* instance_trfs = (Matrix*)std::aligned_alloc(64, meshcount * sizeof(Matrix));
    for(size_t i = 0;i < meshcount;i++){
        double arg = (double(i) / double(meshcount)) * 4 * M_PI;
        instance_trfs[i] = MatrixMultiply(MatrixScale(0.05,0.05,0.05), MatrixTranslate(std::cos(arg) * std::sqrt(arg),-3 + arg/2.0,std::sin(arg) * std::sqrt(arg)));
        //instance_trfs[i] = MatrixMultiply(MatrixRotateX(1.0), MatrixMultiply(MatrixScale(0.2,0.2,0.2), MatrixTranslate(-3,-2,0)));
    }
    //reDrawMeshInstanced(human.meshes[0], instance_trfs, meshcount / 2);
    //reDrawMeshInstanced(human.meshes[0], instance_trfs+5, meshcount / 2);
    for(size_t i = 0;i < meshcount;i++){
        //reDrawMeshTransformed(human.meshes[0], instance_trfs[i]);
    }
    //reDrawSphere(Vector3{0,1000,0}, 3.4f);
    
    reDrawTriangle(Vector3{0,0,0}, Vector3{10,0,0}, Vector3{0,10,10});
    //reDrawTriangle(Vector3{0,0,0}, Vector3{0,0,10}, Vector3{10,0,0});
    
    
    //reDrawSphere(Vector3{2, 0.6, 2}, 2.0f);
    //for(int i = 0;i < cottage.meshCount;i++)
    //    reDrawMeshTransformed(cottage.meshes[i], MatrixTranslate(18.569975, -0.665788, 0.528871));
    //reDrawMeshTransformed(mod.meshes[0], MatrixTranslate(10, 10, 2));
    //reDrawMeshTransformed(mod.meshes[0], MatrixTranslate(10, 10, 2));
    //reDrawMeshTransformed(mod.meshes[0], MatrixTranslate(10, 0, 2));
    //auto vec = reDrawMeshInstanced(human.meshes[0], instance_trfs, meshcount/2);
    //auto vec2 = reDrawMeshInstanced(human.meshes[0], instance_trfs + meshcount/2, meshcount/2);
    //for(auto i : vec2){
    //    std::cout << i << "\n";
    //}
    //return 0;
    //reDrawMeshTransformed(cube.meshes[0], MatrixTranslate(0,0,0));
    //reDrawMeshTransformed(mod.meshes[0], trfs[0]);
    rtcCommitScene(gScene);
    //{
    //    Ray r{.position = Vector3{-3,-3,-3}, .direction = Vector3Normalize(Vector3{1,1.5,1})};
    //    auto x = Raycast(r);
    //    std::cout << x.tfar << "\n";
    //    return 0;
    //}
    Camera3D cam;
    cam.fovy = 45.0f;
    cam.projection = CAMERA_PERSPECTIVE;
    cam.position = Vector3{20, 6, 20};
    cam.target = Vector3{1e-3, 1e-3, 0};
    cam.up = Vector3{0, 1, 0};
    Image rtImage = GenImageColor(GetScreenWidth() / 2, GetScreenHeight(), BLACK);
    Texture tex = LoadTextureFromImage(rtImage);
    int w = GetScreenWidth();
    int wh = GetScreenWidth() / 2;
    int h = GetScreenHeight();
    bool should_raytrace = true;
    bool should_rasterize = true;
    RenderTexture2D rtex = LoadRenderTexture(wh, h);
    Material default_mat = LoadMaterialDefault();
    constexpr int raytrace_supersampling = 1;
    while(!WindowShouldClose()){
        
        ClearBackground(Color { 50, 50, 70, 255 });
        
        if(should_raytrace){
            screenloop([&](int i, int j){
            Vector3 normalized_rgb{0,0,0};
            for(int ssi = 0;ssi < raytrace_supersampling;ssi++){
                for(int ssj = 0;ssj < raytrace_supersampling;ssj++){
                    Ray r = GetScreenToWorldRayEx(Vector2{float(j) + (float(ssj) / raytrace_supersampling), float(i) + (float(ssi) / raytrace_supersampling)}, cam, wh, h);
                    auto hit = Raycast(r);
                    hit.normal = Vector3Normalize(hit.normal);
                    if(!std::isinf(hit.tfar)){
                        normalized_rgb.x += std::clamp(-1.0f * hit.normal.y, 0.0f, 1.0f) * (1.0f / raytrace_supersampling / raytrace_supersampling);
                        normalized_rgb.y += std::clamp(-1.0f * hit.normal.y, 0.0f, 1.0f) * (1.0f / raytrace_supersampling / raytrace_supersampling);
                        normalized_rgb.z += std::clamp(-1.0f * hit.normal.y, 0.0f, 1.0f) * (1.0f / raytrace_supersampling / raytrace_supersampling);
                    }
                }
            }
            ((Color*)rtImage.data)[(i * wh + j)] = Color{
                (unsigned char)(std::min(std::max(normalized_rgb.x,0.0f), 1.0f) * 255),
                (unsigned char)(std::min(std::max(normalized_rgb.y,0.0f), 1.0f) * 255),
                (unsigned char)(std::min(std::max(normalized_rgb.z,0.0f), 1.0f) * 255), 
            255};}, wh, h);
            //{
            //    Ray r{.position = Vector3{-3,-3,-3}, .direction = Vector3Normalize(Vector3{1,1.5,1})};
            //    auto x = Raycast(r);
            //    std::cout << x.tfar << "\n";
            //    return 0;
            //}
            UpdateTexture(tex, rtImage.data);
            DrawTexture(tex,0,0,WHITE);
        }
        
        if(IsKeyPressed(KEY_D)){
            should_rasterize = !should_rasterize;
        }
        if(IsKeyPressed(KEY_R)){
            should_raytrace = !should_raytrace;
        }
        if(should_rasterize){
            BeginTextureMode(rtex);
            ClearBackground(BLACK);
            BeginMode3D(cam);
            for(size_t i = 0;i < meshcount;i++){
                DrawMesh(human.meshes[0], default_mat, instance_trfs[i]);
            }
            //DrawMeshInstanced(human.meshes[0], default_mat, instance_trfs, meshcount);
            EndMode3D();
            EndTextureMode();
            DrawTexturePro(rtex.texture,Rectangle{0,0,(float)wh, (float)-h}, Rectangle{(float)wh, 0, (float)wh, (float)h}, Vector2{0,0},0.0f,WHITE);
        }
        rlSetLineWidth(5.0f);
        DrawLine(wh,0,wh,h,BLACK);
        DrawText(TextFormat("FPS: %d", GetFPS()), 0, 0, 40, GREEN);
        DrawText("Raytraced [Toggle R]", wh / 2 - 100, 0, 40, should_raytrace ? GREEN : RED);
        DrawText("Default [Toggle D]", wh * 3 / 2 - 100, 0, 40, should_rasterize ? GREEN : RED);
        EndDrawing();
    }
}