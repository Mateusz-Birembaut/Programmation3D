#ifndef PHOTONS_H
#define PHOTONS_H

#include "Vec3.h"

struct Photon {
    Vec3 position;  // Position of the photon
    Vec3 direction; // Direction from which the photon came
    Vec3 power;     // Power (color) of the photon

    Photon(const Vec3& pos, const Vec3& dir, const Vec3& pow)
        : position(pos), direction(dir), power(pow) {}

    Photon() {}
};

#endif // PHOTONS_H

/*     std::vector<Photon> photons;
    std::cout << "photon size : " << photons.size() << std::endl;
    std::cout << "Photon map generation started" << std::endl;

    scenes[selected_scene].photonMap(photons); 
    
        Photon photonMap(std::vector<Photon> & photons) {
        for (unsigned int i = 0; i < lights.size(); i++) {
            Light light = lights[i];
            if (light.type == LightType_Spherical) {
                for (int j = 0; j < 1; j++) {
                    Vec3 lightPos = light.pos + Vec3(dist(rng), dist(rng), dist(rng)) * light.radius;
                    Vec3 direction = Vec3(dist(rng), dist(rng), dist(rng));
                    direction.normalize();
                    Ray ray = Ray(lightPos, direction);
                    photonMapRecursive(photons, ray, 3, light);
                }
            }
        }
    }

    void photonMapRecursive(std::vector<Photon> &photons, Ray const &ray, int depth, Light& light) {
        if (depth <= 0) return;
        RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        if (raySceneIntersection.intersectionExists) {
            // Créer un photon à l'endroit de l'intersection
            int index = raySceneIntersection.objectIndex;
            Material mat;
            Vec3 P, N, L, V;
            float u = 0.0f, v = 0.0f;

            if (raySceneIntersection.typeOfIntersectedObject == 0) { // si intersecte une sphere
                P = raySceneIntersection.raySphereIntersection.intersection; 
                N = raySceneIntersection.raySphereIntersection.normal; 

                mat = spheres[index].material;
                u = raySceneIntersection.raySphereIntersection.theta ;
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
            Photon photon;
            photon.position = P;
            photon.direction = ray.direction();
            photon.power = Vec3(0, 1, 0);
            std::cout << "Photon position : " << photon.position << std::endl;
            std::cout << "Photon direction : " << photon.direction << std::endl;
            std::cout << "Photon power : " << photon.power << std::endl;


            photons.push_back(photon);

            Ray reflectedRay =  Ray(P ,reflect(N, ray.direction())); 
            Ray refractedRay = Ray(P ,refract(ray.direction(), N, mat.index_medium)); 

            // Tracer les photons réfléchis et réfractés
            photonMapRecursive(photons, reflectedRay, depth - 1, light);
            photonMapRecursive(photons, refractedRay, depth - 1, light);
        }
    }
    
    */