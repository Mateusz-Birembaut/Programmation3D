
#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H


#include <cmath>
#include <iostream>
#include "Vec3.h"
#include "Line.h"
#include "Ray.h"
#include "Triangle.h"
#include <cfloat>
#include <vector>
#include "Photon.h"


class BoundingBox {
private:

public:
    Vec3 min , max;
    BoundingBox() {}
    BoundingBox( Vec3 const & b , Vec3 const & t )  {
        min = b;
        max = t;
    };

    static BoundingBox photonMapBoundingBox(const std::vector<Photon> & photons) {
        BoundingBox box;
        box.min = Vec3(FLT_MAX, FLT_MAX, FLT_MAX); 
        box.max = Vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX); 

        for ( Photon const & photon : photons ) {
            box.min[0] = std::min(box.min[0], photon.position[0]);
            box.min[1] = std::min(box.min[1], photon.position[1]);
            box.min[2] = std::min(box.min[2], photon.position[2]);

            box.max[0] = std::max(box.max[0], photon.position[0]);
            box.max[1] = std::max(box.max[1], photon.position[1]);
            box.max[2] = std::max(box.max[2], photon.position[2]);
        }

        return box;
    }

    static BoundingBox createSphereBox(const Vec3& center, float radius) {
        BoundingBox box;
        box.min = center - Vec3(radius, radius, radius);
        box.max = center + Vec3(radius, radius, radius);
        return box;
    }

    bool overlaps(const BoundingBox& other) const {
        return (min[0] <= other.max[0] && max[0] >= other.min[0]) &&
               (min[1] <= other.max[1] && max[1] >= other.min[1]) &&
               (min[2] <= other.max[2] && max[2] >= other.min[2]);
    }

    static BoundingBox meshBoundingBox(const std::vector<Triangle>& triangles) {
        BoundingBox box;
        box.min = Vec3(FLT_MAX, FLT_MAX, FLT_MAX); 
        box.max = Vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX); 

        for (const Triangle& triangle : triangles) {
            for (const Vec3& vertex : triangle.getVertices()) { 
                box.min[0] = std::min(box.min[0] , vertex[0] );
                box.min[1] = std::min(box.min[1], vertex[1]);
                box.min[2] = std::min(box.min[2], vertex[2]);

                box.max[0] = std::max(box.max[0], vertex[0]);
                box.max[1] = std::max(box.max[1], vertex[1]);
                box.max[2] = std::max(box.max[2], vertex[2]);
            }
        }

        return box;
    }

    static BoundingBox triangleBoundingBox(const Triangle & triangle) {
        BoundingBox box;
        box.min = Vec3(FLT_MAX, FLT_MAX, FLT_MAX); 
        box.max = Vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX); 
        for (const Vec3& vertex : triangle.getVertices()) { 
            box.min[0] = std::min(box.min[0] , vertex[0] );
            box.min[1] = std::min(box.min[1], vertex[1]);
            box.min[2] = std::min(box.min[2], vertex[2]);

            box.max[0] = std::max(box.max[0], vertex[0]);
            box.max[1] = std::max(box.max[1], vertex[1]);
            box.max[2] = std::max(box.max[2], vertex[2]);
        
        }
        return box;
    }

    std::pair<float, float> intersect(Ray const & ray) const { 

        for(int i = 0; i < 3; i++){ // test si le rayon est parallele a l'axe et n'est pas dans la boite
            if (ray.direction()[i] == 0 && (ray.origin()[i] < min[i] || ray.origin()[i] > max[i])){
                return {INFINITY, INFINITY};
            }
        }

        float rd = 1.0f / ray.direction()[0];

        float tStart = (min[0] - ray.origin()[0]) * rd;
        float tEnd = (max[0] - ray.origin()[0]) * rd;

        if (tStart > tEnd){
            std::swap(tStart, tEnd);
        }

        for(int i = 1; i < 3; i++){
            float rd = 1.0f / ray.direction()[i];

            float t1 = (min[i] - ray.origin()[i]) * rd;
            float t2 = (max[i] - ray.origin()[i]) * rd;

            if (t1 > t2){
                std::swap(t1, t2);
            }

            if (t1 > tStart){
                tStart = t1;
            }
            if (t2 < tEnd){
                tEnd = t2;
            }
        }

        if (tStart > tEnd || tEnd < 0.0001){
            return {INFINITY, INFINITY};
        }

        return {tStart, tEnd};
    }

    static float calculateSurfaceArea(const BoundingBox& box) {
        Vec3 dimensions = box.max - box.min;
        return 2 * (dimensions[0] * dimensions[1] + dimensions[1] * dimensions[2] + dimensions[2] * dimensions[0]);
    }
        

};
#endif 