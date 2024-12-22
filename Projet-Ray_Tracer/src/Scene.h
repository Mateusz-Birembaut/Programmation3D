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


#include <GL/glut.h>


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

    std::vector<KdTreeNode*> kdTrees;

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
    }


    RaySceneIntersection computeIntersection(Ray const & ray) { // recupère l'intersection la plus proche
        RaySceneIntersection result;
        result.t = FLT_MAX; 
        result.intersectionExists = false; 

        for( int It = 0 ; It < spheres.size() ; ++It ) {
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
        for (int i = 0; i < squares.size(); i++){
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

/*      for (int i = 0; i < meshes.size(); i++){
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

       for (int i = 0; i < kdTrees.size(); i++){
            KdTreeNode* kdTree = kdTrees[i];
            RayTriangleIntersection resultTriangleTemp = kdTree->traverse(ray, result.rayMeshIntersection);
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
            
        } 
        
        return result;
    }

    bool intersectObject(Ray const & ray, Vec3 const & lightPos) { // retourne true si un objet est entre le point d'intersection et la lumière et false sinon
        float distToLight = (lightPos - ray.origin()).length();
  
        for( int It = 0 ; It < spheres.size() ; ++It ) {
            Sphere const & sphere = spheres[It];
            RaySphereIntersection resultSphereTemp = sphere.intersect(ray);
            if (resultSphereTemp.intersectionExists ) {
                if (resultSphereTemp.t > 0.001 && resultSphereTemp.t < distToLight) {   
                    if (sphere.material.type == Material_Glass) { // si on intersecte une sphere en verre, on ne le considère pas comme un obstacle
                        continue;
                    }                 
                    return true;
                }
            }
            
        }
        for (int i = 0; i < squares.size(); i++){
            Square const & square = squares[i];
            RaySquareIntersection resultSquareTemp = square.intersect(ray);

            if (resultSquareTemp.intersectionExists){
                if (resultSquareTemp.t > 0.001 && resultSquareTemp.t < distToLight) {
                    return true;
                }
            }
        }

        for (int i = 0; i < meshes.size(); i++){
            Mesh const & m = meshes[i];
            RayTriangleIntersection resultTriangleTemp = m.intersect(ray);

            if (resultTriangleTemp.intersectionExists){
                if (resultTriangleTemp.t > 0.001 && resultTriangleTemp.t < distToLight ) {
                    return true;
                }
            }
        }
        for (int i = 0; i < kdTrees.size(); i++){
            KdTreeNode* kdTree = kdTrees[i];
            RayTriangleIntersection resultTriangleTemp = kdTree->traverse(ray, resultTriangleTemp);
            if (resultTriangleTemp.intersectionExists){
                if (resultTriangleTemp.t > 0.001 && resultTriangleTemp.t < distToLight ) {
                    return true;
                }
            }
            
        } 
        
        return false;
    }


    Vec3 rayTraceRecursive( Ray ray , int NRemainingBounces) {     
        RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        Vec3 color;
        if (raySceneIntersection.intersectionExists) {  
            int index = raySceneIntersection.objectIndex;
            float shininess_obj, transparency_obj, index_medium_obj;
            MaterialType material_type_obj;
            Vec3 P, N, L, V, ambient_obj, diffuse_obj, specular_obj;

            if (raySceneIntersection.typeOfIntersectedObject == 0) { // si intersecte une sphere
                P = raySceneIntersection.raySphereIntersection.intersection; 
                N = raySceneIntersection.raySphereIntersection.normal; 

                extractMaterialProperties(index, spheres[index].material, ambient_obj, diffuse_obj, specular_obj, shininess_obj, material_type_obj, transparency_obj, index_medium_obj);
            } else if (raySceneIntersection.typeOfIntersectedObject == 1) { // si intersecte un carré
                P = raySceneIntersection.raySquareIntersection.intersection; 
                N = raySceneIntersection.raySquareIntersection.normal; 

                extractMaterialProperties(index, squares[index].material, ambient_obj, diffuse_obj, specular_obj, shininess_obj, material_type_obj, transparency_obj, index_medium_obj);
            } else if (raySceneIntersection.typeOfIntersectedObject == 2) { // si intersecte un mesh
                P = raySceneIntersection.rayMeshIntersection.intersection; 
                N = raySceneIntersection.rayMeshIntersection.normal; 
                
                extractMaterialProperties(index, meshes[index].material, ambient_obj, diffuse_obj, specular_obj, shininess_obj, material_type_obj, transparency_obj, index_medium_obj);
            }

            P = P + 0.00001 * N; 

            V = ray.origin() - P; // le vecteur du point d'intersection vers l'origine du rayon
            V.normalize();


            if(material_type_obj != Material_Mirror && material_type_obj != Material_Glass){
            
                for(Light& light : lights) {
                    L = light.pos - P; // le vecteur du point d'intersection vers la lumière
                    L.normalize();

                    Vec3 color_current = shade(light, L, V, N, ray.origin(), shininess_obj, diffuse_obj, specular_obj, ambient_obj);
                    float unblocked_rays_percentage = traceShadowRays(light, P, N, 5);
                    color_current *= unblocked_rays_percentage;

                    color += color_current;
                }

                color /= lights.size();
            }
            
            if( NRemainingBounces > 0 ){
                if (material_type_obj == Material_Mirror ) {
                    Ray reflectedRay =  Ray(P ,reflect(N, ray.direction())); 
                    color = rayTraceRecursive(reflectedRay, NRemainingBounces - 1);
                }else if ( material_type_obj == Material_Glass ) {
                    Ray refractedRay = Ray(P ,refract(ray.direction(), N, index_medium_obj)); 
                    color = rayTraceRecursive(refractedRay, NRemainingBounces - 1);
                }           
            }
            
        } else {
            return Vec3(0.53, 0.81, 0.92);
        }
        return color;
    }

    Vec3 reflect(Vec3 N, Vec3 I) {
        float cosI = -1 * Vec3::dot(N, I);
        Vec3 reflectedDirection = (I + 2 * cosI * N);
        reflectedDirection.normalize();
        return reflectedDirection;
    }

    Vec3 refract(const Vec3 &I, const Vec3 &N,  float eta_t,  float eta_i=1.f) {
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

    Vec3 shade(Light light, Vec3 L, Vec3 V, Vec3 N, Vec3 ray_origin, float shininess_obj , Vec3 diffuse_obj, Vec3 specular_obj, Vec3 ambient_obj) {
        Vec3 R = 2 * Vec3::dot(N, L) * (N) - L;
        R.normalize();

        Vec3 ambient = ambient_obj;
        Vec3 diffuse = diffuse_obj * std::max(0.0f, Vec3::dot(L, N));
        Vec3 specular = specular_obj * std::pow(std::max(0.0f, Vec3::dot(R, V)), shininess_obj);

        return  (ambient + diffuse + specular) * light.material;
    }
    
    float traceShadowRays (Light& light, Vec3 P, Vec3 N, float nbShadowRays) {
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
    

    Vec3 randomPointOnLight(float radius, const Vec3& center) {
        float theta = ((float)rand() / RAND_MAX) * 2 * M_PI;
        float phi = acos(1 - 2 * ((float)rand() / RAND_MAX));
        float x = radius * sin(phi) * cos(theta);
        float y = radius * sin(phi) * sin(theta);
        float z = radius * cos(phi);
        return Vec3(x, y, z) + center; 
    }

    void extractMaterialProperties(int index, Material& material, Vec3& ambient_obj, Vec3& diffuse_obj, Vec3& specular_obj, float& shininess_obj, MaterialType& material_type_obj, float& transparency_obj, float& index_medium_obj) {
        ambient_obj = material.ambient_material;
        diffuse_obj = material.diffuse_material;
        specular_obj = material.specular_material;
        shininess_obj = material.shininess;
        transparency_obj = material.transparency;
        index_medium_obj = material.index_medium;
        material_type_obj = material.type;
    }

    Vec3 rayTrace( Ray const & rayStart ) {
        Vec3 color = rayTraceRecursive(rayStart, 3);
        return color;
    }

    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

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
            s.m_radius = 0.5f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
    }

    void setup_2_spheres() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

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
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.8,0.8,0.8 );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
    }

    void setup_2_planes(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

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
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.03;
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
        }

    }

    void setup_plan_2_spheres(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

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


        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 200., 200.);
            s.translate(Vec3(-50., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.13,1,0.13 );
            s.material.specular_material = Vec3( 1.,1.,1.);
            s.material.shininess = 16;
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

        { // Mesh
            meshes.resize( meshes.size() + 1 );
            Mesh & m = meshes[meshes.size() - 1];
            //m.loadOFF("img/mesh/tetrahedron.off");
            m.loadOFF("img/mesh/nefertiti.off");
            m.translate(Vec3(0., 0., -2.));
            m.build_arrays();
            m.material.diffuse_material = Vec3(0., 0., 1.);
            m.material.specular_material = Vec3(1.0, 1.0, 1.0);
            m.material.shininess = 32;

            // Extraire les triangles du maillage
            std::vector<Triangle> triangles;
            for (const auto& meshTriangle : m.triangles) {
                Vec3 v0 = m.vertices[meshTriangle[0]].position;
                Vec3 v1 = m.vertices[meshTriangle[1]].position;
                Vec3 v2 = m.vertices[meshTriangle[2]].position;
                triangles.emplace_back(v0, v1, v2);
            }

            KdTreeNode* kdTree = new KdTreeNode();
            //kdTree->KdTreeBuild(triangles, 10);
            //kdTrees.push_back(kdTree);
            //kdTree->testKdTreeBuild();
            //kdTree->testTraverse();
            kdTree->testKdTreeBuildComplex();
        }

    }

};

#endif
