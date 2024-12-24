
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


class BoundingBox {
private:

public:
    Vec3 min , max;
    BoundingBox() {}
    BoundingBox( Vec3 const & b , Vec3 const & t )  {
        min = b;
        max = t;
    };


    BoundingBox meshBoundingBox(const std::vector<Triangle>& triangles) {
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

    BoundingBox triangleBoundingBox(const Triangle & triangle) {
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
        

};
#endif 