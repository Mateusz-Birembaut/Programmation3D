#ifndef Sphere_H
#define Sphere_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySphereIntersection{
    bool intersectionExists;
    float t;
    float theta,phi;
    Vec3 intersection;
    Vec3 secondintersection;
    Vec3 normal;
};


static
Vec3 SphericalCoordinatesToEuclidean( float theta , float phi ) {
    return Vec3( cos(theta) * cos(phi) , sin(theta) * cos(phi) , sin(phi) );
}

static
Vec3 SphericalCoordinatesToEuclidean( Vec3 ThetaPhiR ) {
    return ThetaPhiR[2] * Vec3( cos(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[1]) );
}

static
Vec3 EuclideanCoordinatesToSpherical( Vec3 xyz ) {
    float R = xyz.length();
    float phi = asin( xyz[2] / R );
    float theta = atan2( xyz[1] , xyz[0] );
    return Vec3( theta , phi , R );
}



class Sphere : public Mesh {
public:
    Vec3 m_center;
    float m_radius;

    Sphere() : Mesh() {}
    Sphere(Vec3 c , float r) : Mesh() , m_center(c) , m_radius(r) {}

    void build_arrays(){
        unsigned int nTheta = 20 , nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi );
        normalsArray.resize(3 * nTheta * nPhi );
        uvs_array.resize(2 * nTheta * nPhi );
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta ; ++thetaIt ) {
            float u = (float)(thetaIt) / (float)(nTheta-1);
            float theta = u * 2 * M_PI;
            for( unsigned int phiIt = 0 ; phiIt < nPhi ; ++phiIt ) {
                unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                float v = (float)(phiIt) / (float)(nPhi-1);
                float phi = - M_PI/2.0 + v * M_PI;
                Vec3 xyz = SphericalCoordinatesToEuclidean( theta , phi );
                positions_array[ 3 * vertexIndex + 0 ] = m_center[0] + m_radius * xyz[0];
                positions_array[ 3 * vertexIndex + 1 ] = m_center[1] + m_radius * xyz[1];
                positions_array[ 3 * vertexIndex + 2 ] = m_center[2] + m_radius * xyz[2];
                normalsArray[ 3 * vertexIndex + 0 ] = xyz[0];
                normalsArray[ 3 * vertexIndex + 1 ] = xyz[1];
                normalsArray[ 3 * vertexIndex + 2 ] = xyz[2];
                uvs_array[ 2 * vertexIndex + 0 ] = u;
                uvs_array[ 2 * vertexIndex + 1 ] = v;
            }
        }
        triangles_array.clear();
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta - 1 ; ++thetaIt ) {
            for( unsigned int phiIt = 0 ; phiIt < nPhi - 1 ; ++phiIt ) {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt+1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt+1) * nTheta;
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuV );
            }
        }
    }


    Vec3 getSphereUV(const Vec3& p) const{
        float u = 0.5f + atan2(p[2], p[0]) / (2.0f * M_PI);
        float v = 0.5f - asin(p[1]) / M_PI;
        return Vec3(u, v, 0);
    }


    RaySphereIntersection intersect(const Ray &ray) const {
        RaySphereIntersection intersection;
        intersection.intersectionExists = false;

        Vec3 ray_origin = ray.origin();
        Vec3 ray_direction = ray.direction();

        Vec3 origin_center =  ray_origin - this->m_center;
        float origin_center_length = origin_center.length();

        float b = 2.0f * (Vec3::dot(ray_direction , origin_center));
        float c = origin_center_length * origin_center_length  - this->m_radius * this->m_radius;

        float discriminant = b*b - 4 * c ;

        if (discriminant < 0) {
            return intersection;
        }

        if (discriminant == 0) {
            float t = -b / 2.0;
            if (t > 0){
                intersection.intersection = ray_origin + t * ray_direction;
                intersection.intersectionExists = true;
                intersection.t = t;
                Vec3 normal = intersection.intersection - this->m_center;
                normal.normalize();
                intersection.normal = normal;
                return intersection;
            }else {
                return intersection;
            }

        }

        float sqrtDiscriminant = sqrt(discriminant);

        float t1 = (-b + sqrtDiscriminant) / 2.0;
        float t2 = (-b - sqrtDiscriminant) / 2.0;

        float min_positif = getSmallestT(t1, t2);

        if (min_positif == -1){
            return intersection;
        }
        
        Vec3 const intersection_point = ray_origin + min_positif * ray_direction;

        intersection.intersectionExists = true;
        intersection.t = min_positif;
        intersection.intersection = intersection_point;

        Vec3 normal = (intersection_point - this->m_center);
        normal.normalize();
        intersection.normal = normal;

        Vec3 uv = getSphereUV(normal);

        intersection.theta = uv[0];
        intersection.phi = uv[1];

        return intersection;
    }

    float getSmallestT(float t1, float t2) const {
        if (t1 < 0 && t2 < 0) {
            return -1;
        }
        if (t1 < 0) {
            return t2;
        }
        if (t2 < 0) {
            return t1;
        }
        return std::min(t1, t2);

    }


};
#endif
