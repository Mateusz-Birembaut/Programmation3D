#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"
#include "Ray.h"
#include "Sphere.h"
#include "Square.h"
#include "Vec3.h"
#include "BoundingBox.h"
#include "KdTree.h"
#include "Photon.h"


#include <GL/glut.h>

#include <random>

std::random_device rand_dev;
std::mt19937 rng(rand_dev());
static std::uniform_real_distribution<float> dist(0.0f, 1.0f);
static std::uniform_real_distribution<float> dist11(-1.0f, 1.0f);
static std::uniform_real_distribution<float> dist05(-0.5f, 0.5f);

enum LightType {
    LightType_Spherical,
    LightType_Quad
};


struct Light {
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Mesh quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}

};

struct RaySceneIntersection{
    bool intersectionExists;
    unsigned int typeOfIntersectedObject;
    unsigned int objectIndex;
    float t;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;
    RaySceneIntersection() : intersectionExists(false) , t(FLT_MAX) {}
};



class Scene {

    std::vector< Mesh > meshes;
    std::vector< Sphere > spheres;
    std::vector< Square > squares;
    std::vector< Light > lights;

    std::vector<KdTree> kdTrees;
    std::vector<ppmLoader::ImageRGB> textures;
    std::vector<ppmLoader::ImageRGB> normalMaps;

public:

    Scene() {
    }

    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
         for( unsigned int It = 0 ; It < meshes.size() ; ++It ) {
            Mesh const & mesh = meshes[It];
            mesh.draw();
        } 
        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) {
            Sphere const & sphere = spheres[It];
            sphere.draw();
        }
        for( unsigned int It = 0 ; It < squares.size() ; ++It ) {
            Square const & square = squares[It];
            square.draw();
        }
        for( unsigned int It = 0 ; It < kdTrees.size() ; ++It ) {
            //KdTree kdTree = kdTrees[It];
            //kdTree.drawBoundingBoxesHelper(kdTree.root, false);
        }
    }


    RaySceneIntersection computeIntersection(Ray const & ray) { // recupère l'intersection la plus proche
        RaySceneIntersection result;
        result.t = FLT_MAX; 
        result.intersectionExists = false; 

        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) {
            Sphere const & sphere = spheres[It];
            RaySphereIntersection resultSphereTemp = sphere.intersect(ray);
            if (resultSphereTemp.intersectionExists ) {
                if (resultSphereTemp.t > 0.001 && resultSphereTemp.t < result.t) {                  
                    result.t = resultSphereTemp.t;
                    result.intersectionExists = true;
                    result.objectIndex = It;
                    result.raySphereIntersection = resultSphereTemp;
                    result.typeOfIntersectedObject = 0;
                }
            }
            
        }
        for (unsigned int i = 0; i < squares.size(); i++){
            Square const & square = squares[i];
            RaySquareIntersection resultSquareTemp = square.intersect(ray);

            if (resultSquareTemp.intersectionExists){
                if (resultSquareTemp.t > 0.001 && resultSquareTemp.t < result.t ) {
                    result.intersectionExists = true;
                    result.t = resultSquareTemp.t;
                    result.objectIndex = i;
                    result.raySquareIntersection = resultSquareTemp;
                    result.typeOfIntersectedObject = 1;
                    result.raySquareIntersection.normal = resultSquareTemp.normal;
                }
            }
        }

/*      for (unsigned int i = 0; i < meshes.size(); i++){
            Mesh const & m = meshes[i];

            RayTriangleIntersection resultTriangleTemp = m.intersect(ray);

            if (resultTriangleTemp.intersectionExists){
                if (resultTriangleTemp.t > 0.001 && resultTriangleTemp.t < result.t ) {
                    result.intersectionExists = true;
                    result.t = resultTriangleTemp.t;
                    result.objectIndex = i;
                    result.rayMeshIntersection = resultTriangleTemp;
                    result.typeOfIntersectedObject = 2;
                    result.raySquareIntersection.normal = resultTriangleTemp.normal;
                }
            }
        } */

        for( unsigned int i = 0; i < kdTrees.size(); i++){
            KdTree kdTree = kdTrees[i];
            std::pair<float, float> interval = kdTree.root->node_box.intersect(ray);

            if (interval.first == INFINITY && interval.second == INFINITY){// si le rayon ne touche pas la boite
                continue;
            }

            RayTriangleIntersection resultTriangleTemp = kdTree.traverse(ray, kdTree.root, interval.first, interval.second);
            if (resultTriangleTemp.intersectionExists){
                if (resultTriangleTemp.t > 0.001 && resultTriangleTemp.t < result.t ) {
                    result.intersectionExists = true;
                    result.t = resultTriangleTemp.t;
                    result.objectIndex = i;
                    result.rayMeshIntersection = resultTriangleTemp;
                    result.typeOfIntersectedObject = 2;
                    result.rayMeshIntersection.normal = resultTriangleTemp.normal;
                }
            }
        }
        
        return result;
    }

    bool intersectObject(Ray const & ray, Vec3 const & lightPos) { // retourne true si un objet est entre le point d'intersection et la lumière et false sinon
        float distToLight = (lightPos - ray.origin()).length();
  
        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) {
            Sphere const & sphere = spheres[It];
            RaySphereIntersection resultSphereTemp = sphere.intersect(ray);
            if (resultSphereTemp.intersectionExists ) {
                if (resultSphereTemp.t > 0.001 && resultSphereTemp.t < distToLight) {   
                    /*if (sphere.material.type == Material_Glass) { // si on intersecte une sphere en verre, on ne le considère pas comme un obstacle
                        continue;
                    }*/               
                    return true;
                }
            }
            
        }
        for (unsigned int i = 0; i < squares.size(); i++){
            Square const & square = squares[i];
            RaySquareIntersection resultSquareTemp = square.intersect(ray);

            if (resultSquareTemp.intersectionExists){
                if (resultSquareTemp.t > 0.001 && resultSquareTemp.t < distToLight) {
                    return true;
                }
            }
        }

/*         for (int i = 0; i < meshes.size(); i++){
            Mesh const & m = meshes[i];
            RayTriangleIntersection resultTriangleTemp = m.intersect(ray);

            if (resultTriangleTemp.intersectionExists){
                if (resultTriangleTemp.t > 0.001 && resultTriangleTemp.t < distToLight ) {
                    return true;
                }
            }
        }  */

       for (unsigned int i = 0; i < kdTrees.size(); i++){
            KdTree kdTree = kdTrees[i];
            std::pair<float, float> interval = kdTree.root->node_box.intersect(ray);

            if (interval.first == INFINITY && interval.second == INFINITY){// si le rayon ne touche pas la boite
                continue;
            }

            RayTriangleIntersection resultTriangleTemp = kdTree.traverse(ray, kdTree.root, interval.first, interval.second);
            if (resultTriangleTemp.intersectionExists){
                if (resultTriangleTemp.t > 0.001 && resultTriangleTemp.t < distToLight ) {
                    return true;
                }
            }
        }

        return false;
    }


    void handleIntersection(const RaySceneIntersection& raySceneIntersection,const Ray & ray , Vec3& P, Vec3& N, Vec3& V, Material& mat, float& u, float& v, int index) {
        if (raySceneIntersection.typeOfIntersectedObject == 0) { // si intersecte une sphere
            P = raySceneIntersection.raySphereIntersection.intersection; 
            N = raySceneIntersection.raySphereIntersection.normal; 

            mat = spheres[index].material;
            u = raySceneIntersection.raySphereIntersection.theta;
            v = raySceneIntersection.raySphereIntersection.phi;

        } else if (raySceneIntersection.typeOfIntersectedObject == 1) { // si intersecte un carré
            P = raySceneIntersection.raySquareIntersection.intersection; 
            N = raySceneIntersection.raySquareIntersection.normal; 

            mat = squares[index].material;
            u = raySceneIntersection.raySquareIntersection.u;
            v = raySceneIntersection.raySquareIntersection.v;
            
        } else if (raySceneIntersection.typeOfIntersectedObject == 2) { // si intersecte un mesh
            P = raySceneIntersection.rayMeshIntersection.intersection; 
            N = raySceneIntersection.rayMeshIntersection.normal; 
            
            mat = meshes[index].material;
        }

        if( mat.texture != nullptr ) {
            mat.diffuse_material = mat.texture->sampleTextureAsVec3(u,v, mat.t_uRepeat, mat.t_vRepeat);
        }   

        if (mat.normalMap != nullptr) {
            updateNormal(N, mat.normalMap, u, v, mat, raySceneIntersection);
        } 

        P = P + 0.00001 * N; // decale le point d'intersection pour eviter de le reintersecter

        V = ray.origin() - P; // le vecteur du point d'intersection vers l'origine du rayon
        V.normalize();

    }

    Vec3 rayTraceRecursive(KdTreePhotonMap & photonMap, Ray ray , int NRemainingBounces) {     
        RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        Vec3 color;
        if (raySceneIntersection.intersectionExists) {  
            int index = raySceneIntersection.objectIndex;
            Material mat;
            Vec3 P, N, L, V;
            float u = 0.0f, v = 0.0f;

            handleIntersection(raySceneIntersection, ray, P, N, V, mat, u, v, index);

            
            if(mat.type != Material_Mirror && mat.type  != Material_Glass){
            
                Vec3 directLight(0,0,0);
                for(Light& light : lights) {
                    L = light.pos - P; // le vecteur du point d'intersection vers la lumière
                    L.normalize();

                    float unblocked_rays_percentage = traceShadowRays(light, P, 5);
                    if (unblocked_rays_percentage == 0){
                        continue;
                    }
                    
                    Vec3 color_current = shade(light, L, V, N, ray.origin(), mat, u , v);
                    directLight += color_current * unblocked_rays_percentage;
                }

                directLight /= lights.size();

                //Vec3 indirectLight = Vec3(0,0,0);
                //indirectLight = estimateRadiance(P, N, photonMap, 0.2f);
                
                //color = indirectLight / 800 ; //directLight + indirectLight; // Scale indirect lighting
            
                color = directLight;
            }
                 
            if( NRemainingBounces > 0 ){
                if (mat.type == Material_Mirror ) {
                    Ray reflectedRay =  Ray(P ,reflect(N, ray.direction())); 
                    color += rayTraceRecursive(photonMap ,reflectedRay, NRemainingBounces - 1);
                }else if ( mat.type == Material_Glass ) {
                    Ray refractedRay = Ray(P ,refract(ray.direction(), N, mat.index_medium)); 
                    color += rayTraceRecursive(photonMap ,refractedRay, NRemainingBounces - 1);
                }           
            }
            
            
        } else {
            return Vec3(0.53, 0.81, 0.92);
        }
        return color;
    }

    Vec3 estimateRadiance(const Vec3& point, const Vec3& normal, KdTreePhotonMap& kdTreePhoton, float radius) {
        Vec3 radiance(0,0,0);
        int photonCount = 0;
        float maxRadius2 = radius * radius;

        std::vector<Photon> nearbyPhotons = kdTreePhoton.findNearestPhotons(point, radius);
        
        for(const Photon& photon : nearbyPhotons) {
            float dist2 = (point - photon.position).length();
            dist2 = dist2 * dist2;
            if(dist2 < maxRadius2) {
                float weight = 1.0f - sqrt(dist2)/radius;
                float cosTheta = Vec3::dot(normal, photon.direction);
                if(cosTheta > 0) {
                    radiance += photon.power * weight * cosTheta;
                    photonCount++;
                }
            }
        }
        
        if(photonCount > 0) {
            float area = M_PI * radius * radius;
            radiance *= (1.0f / area);
        }
        
        return radiance;
    }


    Vec3 transformToWorldSpace(const Vec3& tangentSpaceNormal, const Mat3& TBN) {
        Vec3 worldNormal;
        worldNormal[0] = TBN(0, 0) * tangentSpaceNormal[0] + TBN(0, 1) * tangentSpaceNormal[1] + TBN(0, 2) * tangentSpaceNormal[2];
        worldNormal[1] = TBN(1, 0) * tangentSpaceNormal[0] + TBN(1, 1) * tangentSpaceNormal[1] + TBN(1, 2) * tangentSpaceNormal[2];
        worldNormal[2] = TBN(2, 0) * tangentSpaceNormal[0] + TBN(2, 1) * tangentSpaceNormal[1] + TBN(2, 2) * tangentSpaceNormal[2];
        worldNormal.normalize();
        return worldNormal;
    }

    Mat3 computeTBN(const Vec3& tangent, const Vec3& bitangent, const Vec3& normal) {
        return Mat3(tangent, bitangent, normal);
    }

    float sphereSurfaceArea(float radius) {
        return 4 * M_PI * (radius * radius);
    }

    Vec3 reflect(Vec3& N, const Vec3& I) {
        float cosI = -1 * Vec3::dot(N, I);
        Vec3 reflectedDirection = (I + 2 * cosI * N);
        reflectedDirection.normalize();
        return reflectedDirection;
    }

    Vec3 refract(const Vec3 &I, Vec3 &N,  float eta_t,  float eta_i=1.f) {
        float cosi = -std::max(-1.f, std::min(1.f, Vec3::dot(I, N)));
        Vec3 n = N;
        if (cosi<0){
            n = -1 * N; 
            std::swap(eta_i, eta_t);
        }
        float eta = eta_i / eta_t;
        float k = 1 - eta*eta*(1 - cosi*cosi);
        if (k < 0){
            return reflect(N,I);
        }
        Vec3 unit_vector = I*eta + n*(eta*cosi - sqrtf(k));
        unit_vector.normalize();
        return unit_vector;
    }

    Vec3 shade(Light& light, Vec3& L, Vec3& V, Vec3& N, const Vec3& ray_origin, Material& mat, float& u, float& v) {
        Vec3 R = 2 * Vec3::dot(N, L) * (N) - L;
        R.normalize();

        Vec3 ambient = mat.ambient_material;

        Vec3 diffuse = mat.diffuse_material * std::max(0.0f, Vec3::dot(L, N));

        Vec3 specular = mat.specular_material * std::pow(std::max(0.0f, Vec3::dot(R, V)), mat.shininess);

        return  (ambient + diffuse + specular) * light.material;
    }
    
    float traceShadowRays (Light& light, Vec3& P, float nbShadowRays) {
        int number_shadow_rays = 0;

        for (int i = 0; i < nbShadowRays; i++) {
            Vec3  offseted_light_pos = randomPointOnLight(light.radius, light.pos);
            Vec3 direction = offseted_light_pos - P;
            direction.normalize();

            if (intersectObject(Ray(P, direction), offseted_light_pos)) {
                number_shadow_rays++;
            }
        }
        return 1 - (number_shadow_rays / nbShadowRays); // retourne le pourcentage de rayons non bloqués
    }
    

    Vec3 randomPointOnLight(float& radius, const Vec3& center) {
        float u = dist(rng);
        float v = dist(rng);
        float theta = u * 2.0f * M_PI;
        float phi = acos(1.0f - 2.0f * v);
        float x = radius * sin(phi) * cos(theta);
        float y = radius * sin(phi) * sin(theta);
        float z = radius * cos(phi);
        return Vec3(x, y, z) + center;
    }

    Vec3 rayTrace( KdTreePhotonMap & photonMap,Ray const & rayStart ) {
        Vec3 color = rayTraceRecursive(photonMap ,rayStart, 3);
        return color;
    }


    void photonMap(std::vector<Photon> & photons, unsigned int rayons) {
        for (unsigned int i = 0; i < lights.size(); i++) {
            Light light = lights[i];
            if (light.type == LightType_Spherical) {
                for (unsigned int j = 0; j < rayons; j++) {
                    Vec3 lightPos = randomPointOnLight(light.radius, light.pos);
                    Vec3 direction = Vec3(dist11(rng), dist11(rng), dist11(rng));
                    direction.normalize();
                    Ray ray = Ray(lightPos, direction);
                    photonMapRecursive(photons, ray, 4, 4, light.material, light);
                }
            }
        }
    }

    void photonMapRecursive(std::vector<Photon>& photons, Ray const& ray, int maxDepth, int depth, Vec3& power, Light& light) {
        if (depth <= 0) return;

        RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        if (!raySceneIntersection.intersectionExists) return;

        int index = raySceneIntersection.objectIndex;
        Material mat;
        Vec3 P, N, V;
        float u = 0.0f, v = 0.0f;
        handleIntersection(raySceneIntersection, ray, P, N, V, mat, u, v, index);

        // Store photon on ALL surfaces except perfect mirrors
        Photon photon(P, -1 * ray.direction(), power, N, maxDepth - depth, 2.0f);
        
        // Apply material color to photon power
        photon.power = photon.power * mat.diffuse_material;
        
        if (photon.isSignificant()) {
            photons.push_back(photon);
        }

        // Different bounce behaviors
        if (mat.type == Material_Mirror) {
            Ray reflectedRay = Ray(P, reflect(N, ray.direction()));
            Vec3 newPower = power * 0.9f; // Less energy loss for mirrors
            photonMapRecursive(photons, reflectedRay, maxDepth, depth - 1, newPower, light);
        }
        else if (mat.type == Material_Glass) {
            Vec3 newPower = power * 0.8f;
            Ray reflectedRay = Ray(P, reflect(N, ray.direction()));
            Ray refractedRay = Ray(P, refract(ray.direction(), N, mat.index_medium));
            
            photonMapRecursive(photons, reflectedRay, maxDepth, depth - 1, newPower, light);
            photonMapRecursive(photons, refractedRay, maxDepth, depth - 1, newPower, light);
        }
        else {
            // Diffuse bounce (optional for better GI)
            Vec3 newPower = power * 0.5f;
            Vec3 randomDir = Vec3(dist11(rng), dist11(rng), dist11(rng));
            randomDir.normalize();
            if (Vec3::dot(randomDir, N) < 0) randomDir = -1 * randomDir;
            Ray diffuseRay = Ray(P, randomDir);
            photonMapRecursive(photons, diffuseRay, maxDepth, depth - 1, newPower, light);
        }
    }


    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        kdTrees.clear();
        textures.clear();
        normalMaps.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize( spheres.size() + 2 );
            Sphere & s = spheres[spheres.size() - 2];
            s.m_center = Vec3(1.0,0.0,0.0);
            s.m_radius = 2.f;
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;

            auto texture = new ppmLoader::ImageRGB;
            ppmLoader::load_ppm(*texture, "img/textures/damier.ppm");
            s.material.texture = texture;

            auto normalMap = new ppmLoader::ImageRGB;
            ppmLoader::load_ppm(*normalMap, "img/normalMaps/n1.ppm");
            s.material.normalMap = normalMap;

        }
    }

    void updateNormal(Vec3 & N, ppmLoader::ImageRGB * normalMap, float u, float v, Material & mat, const RaySceneIntersection & raySceneIntersection) {
        Vec3 tangentSpaceNormal = normalMap->sampleNormMapAsVec3(u, v, mat.n_uRepeat, mat.n_vRepeat);
        tangentSpaceNormal = 2 * tangentSpaceNormal - Vec3(1, 1, 1); // Conversion de [0,1] à [-1,1]
        tangentSpaceNormal.normalize();                               
                              
        Vec3 tangent;
        if( raySceneIntersection.typeOfIntersectedObject == 0){
            tangent = Vec3(-sin(v), 0, cos(u));
            tangent.normalize();
        }else if (raySceneIntersection.typeOfIntersectedObject == 1){
            tangent = Vec3(0.0f, 1.0f, 0.0f);
        }

        Vec3 bitangent = Vec3::cross(N, tangent);
        bitangent.normalize();

        // Calcul de la matrice TBN
        Mat3 TBN = computeTBN(tangent, bitangent, N);

        // Transformation vers l'espace monde
        Vec3 worldNormal = transformToWorldSpace(tangentSpaceNormal, TBN); 
        N = worldNormal; // Remplace la normale géométrique par la normale de la normal map
    }

    void setup_2_spheres() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        kdTrees.clear();
        textures.clear();
        normalMaps.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize(spheres.size() + 2); 

            // Première sphère
            Sphere & s1 = spheres[spheres.size() - 1]; 
            s1.m_center = Vec3(1.0, 0.0, 0.0);
            s1.m_radius = 0.5f;
            s1.build_arrays();
            s1.material.type = Material_Mirror;
            s1.material.diffuse_material = Vec3(1., 0., 0.);
            s1.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s1.material.shininess = 20;

            // Deuxième sphère
            Sphere & s2 = spheres[spheres.size() - 2]; 
            s2.m_center = Vec3(-1.0, 0.0, 0.0);
            s2.m_radius = 0.5f;
            s2.build_arrays();
            s2.material.type = Material_Mirror;
            s2.material.diffuse_material = Vec3(0., 1., 0.);
            s2.material.specular_material = Vec3(0.2, 0.2, 0.2);
            s2.material.shininess = 20;

            
        }
    }

    void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        kdTrees.clear();
        textures.clear();
        normalMaps.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        
        {
            textures.resize( textures.size() + 1 );
            ppmLoader::ImageRGB & texture = textures[textures.size() - 1];
            ppmLoader::load_ppm(texture, "img/textures/s6.ppm");
        }

        {
            normalMaps.resize( normalMaps.size() + 1 );
            ppmLoader::ImageRGB & normalMap = normalMaps[normalMaps.size() - 1];
            ppmLoader::load_ppm(normalMap, "img/normalMaps/n4.ppm");
        }

        {
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.8,0.8,0.8 );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;

            s.material.texture = &textures[textures.size() - 1];
            s.material.normalMap = &normalMaps[normalMaps.size() - 1];

        }
    }

    void setup_2_planes(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        kdTrees.clear();
        textures.clear();
        normalMaps.clear();
        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(0,1,1);
            light.isInCamSpace = false;
        }
        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }
    
    }

    void setup_cornell_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        kdTrees.clear();
        textures.clear();
        normalMaps.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 0.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        {
            textures.resize( textures.size() + 1 );
            ppmLoader::ImageRGB & texture = textures[textures.size() - 1];
            ppmLoader::load_ppm(texture, "img/textures/s6.ppm");
        }

        {
            normalMaps.resize( normalMaps.size() + 1 );
            ppmLoader::ImageRGB & normalMap = normalMaps[normalMaps.size() - 1];
            ppmLoader::load_ppm(normalMap, "img/normalMaps/n2.ppm");
        }

        {
            normalMaps.resize( normalMaps.size() + 1 );
            ppmLoader::ImageRGB & normalMap = normalMaps[normalMaps.size() - 1];
            ppmLoader::load_ppm(normalMap, "img/normalMaps/n1.ppm");
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.5,0.5,0.5 );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,1.,1.);
            s.material.shininess = 16;

        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;


        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,0.5 );
            s.material.specular_material = Vec3( 1.,1.,1.);
            s.material.shininess = 16;

            //s.material.texture = &textures[textures.size() - 1];
/*             s.material.normalMap = &normalMaps[normalMaps.size() - 1];
            s.material.n_uRepeat = 1;
            s.material.n_vRepeat = 1; */
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.5,1.0,1.0 );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
        }

        { //Front Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,0.0,1.0 );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
        }

        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            //s.m_center = Vec3(0.0, 0., 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.03;//1.03;

            s.material.normalMap = &normalMaps[normalMaps.size() - 1];
            s.material.n_uRepeat = 3;
            s.material.n_vRepeat = 3;

        }


        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;

/*             s.material.normalMap = &normalMaps[normalMaps.size() - 1];
            s.material.n_uRepeat = 4;
            s.material.n_vRepeat = 4; */
        }
        {
            meshes.resize( meshes.size() + 1 );
            Mesh & m = meshes[meshes.size() - 1];
            //m.loadOFF("img/mesh/turtle.off");
            m.loadOFF("img/mesh/nefertiti.off");
            //m.loadOFF("img/mesh/bimba_3.7Mf.off");
            m.translate(Vec3(-1., 0., -1.));
            m.build_arrays();
            m.material.diffuse_material = Vec3(0., 0., 1.);
            m.material.specular_material = Vec3(1.0, 1.0, 1.0);
            m.material.shininess = 32;

            KdTree kdTree( &m ,8);
            kdTrees.push_back(kdTree);


        }

        {
            meshes.resize( meshes.size() + 1 );
            Mesh & m = meshes[meshes.size() - 1];
            //m.loadOFF("img/mesh/turtle.off");
            m.loadOFF("img/mesh/nefertiti.off");
            //m.loadOFF("img/mesh/bimba_3.7Mf.off");
            m.translate(Vec3(1., 0., -1.));
            m.build_arrays();
            m.material.diffuse_material = Vec3(0., 0., 1.);
            m.material.specular_material = Vec3(1.0, 1.0, 1.0);
            m.material.shininess = 32;

            KdTree kdTree( &m ,8);
            kdTrees.push_back(kdTree);


        }

    }

    void setup_plan_2_spheres(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();
        kdTrees.clear();
        textures.clear();
        normalMaps.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 100.0, 100.0, 10.0 );
            light.radius = 50.0f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        {
            textures.resize( textures.size() + 1 );
            auto & texture = textures[textures.size() - 1];
            ppmLoader::load_ppm(texture, "img/textures/damier.ppm");
        }

        {
            textures.resize( textures.size() + 1 );
            auto & texture = textures[textures.size() - 1];
            ppmLoader::load_ppm(texture, "img/textures/s1.ppm");
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 200., 200.);
            s.translate(Vec3(-50., -50., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.13,1,0.13 );
            s.material.specular_material = Vec3( 1.,1.,1.);
            s.material.shininess = 16;

            s.material.texture = & textures[ textures.size() - 2 ];
            s.material.t_uRepeat = 50;
            s.material.t_vRepeat = 100;
            
        }


        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }


        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.rotate_x(90);
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0.0, -1.25, 2.5);
            s.m_radius = 0.75f;
            s.scale(Vec3(3., 3., 3.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
            

            s.material.texture = & textures[ textures.size() - 1 ];


        }

        { // Mesh
            meshes.resize( meshes.size() + 1 );
            Mesh & m = meshes[meshes.size() - 1];
            //m.loadOFF("img/mesh/turtle.off");
            m.loadOFF("img/mesh/nefertiti.off");
            //m.loadOFF("img/mesh/bimba_3.7Mf.off");
            m.translate(Vec3(0., 0., -2.));
            m.build_arrays();
            m.material.diffuse_material = Vec3(0., 0., 1.);
            m.material.specular_material = Vec3(1.0, 1.0, 1.0);
            m.material.shininess = 32;

            KdTree kdTree( &m ,8);
            kdTrees.push_back(kdTree);

            // stocke tout dans le kdTree
            //peut faire ça pour avoir plus de mémoire
            /*          
            m.triangles.clear();
            m.triangles_array.clear();
            m.vertices.clear(); 
            m.positions_array.clear();
            m.normalsArray.clear();
            m.uvs_array.clear(); 
            */
            

        }

    }

};

#endif
