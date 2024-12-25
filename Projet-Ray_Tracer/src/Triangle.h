#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Vec3.h"
#include "Ray.h"
#include "Plane.h"

struct RayTriangleIntersection{
    bool intersectionExists = false;
    float t;
    float w0,w1,w2;
    unsigned int tIndex;
    Vec3 intersection;
    Vec3 normal;
};

class Triangle {
private:
    Vec3 m_c[3];
    Vec3 m_normal;
    float area;
    float m_D;
public:
    Triangle() {}
    Triangle( Vec3 const & c0 , Vec3 const & c1 , Vec3 const & c2 ) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
    }
    void updateAreaAndNormal() {
        Vec3 nNotNormalized = Vec3::cross( m_c[1] - m_c[0] , m_c[2] - m_c[0] );
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
        m_D = -Vec3::dot(m_normal, m_c[0]);
    }

    float getArea() const { return area; }
    
    void setC0( Vec3 const & c0 ) { m_c[0] = c0; } // remember to update the area and normal afterwards!
    void setC1( Vec3 const & c1 ) { m_c[1] = c1; } // remember to update the area and normal afterwards!
    void setC2( Vec3 const & c2 ) { m_c[2] = c2; } // remember to update the area and normal afterwards!

    std::vector<Vec3> getVertices() const {return {m_c[0], m_c[1], m_c[2]};}

    Vec3 getCenter() const { return (m_c[0] + m_c[1] + m_c[2]) / 3.f; }

    Vec3 const & normal() const { return m_normal; }

    Vec3 projectOnSupportPlane( Vec3 const & p ) const {
        Vec3 result;
        //TODO completer
        return result;
    }

    float squareDistanceToSupportPlane( Vec3 const & p ) const {
        float result;
        //TODO completer
        return result;
    }

    float distanceToSupportPlane( Vec3 const & p ) const { return sqrt( squareDistanceToSupportPlane(p) ); }

    bool isParallelTo( Line const & L ) const {
        if (Vec3::dot(L.direction(), m_normal) >= 0.0001){
            return true;
        }
        return false;
    }
    Vec3 getIntersectionPointWithSupportPlane( Line const & L ) const {
        // you should check first that the line is not parallel to the plane!
        Vec3 result;
        //TODO completer
        return result;
    }
    void computeBarycentricCoordinates( Vec3 const & p , float & u0 , float & u1 , float & u2 ) const {
        //TODO Complete
    }

    RayTriangleIntersection getIntersection( Ray const & ray ) const {
        RayTriangleIntersection result;

        // 1) check that the ray is not parallel to the triangle:
        if(isParallelTo(ray)){
            //std::cout << "le triangle est parallel to the triangle" << std::endl;
            return result;
        }

        Vec3 e1 = m_c[1] - m_c[0];
        Vec3 e2 = m_c[2] - m_c[0];

        // 2) check that the triangle is "in front of" the ray:
        float t = - (Vec3::dot(m_normal, ray.origin()) + m_D) / Vec3::dot(m_normal, ray.direction());
        if (t <= 0.001){
            //std::cout << "le rayon est derriere le rayon" << std::endl;
            return result;
        }

        // 3) check that the intersection point is inside the triangle:
        // CONVENTION: compute u,v such that p = w0*c0 + w1*c1 + w2*c2, check that 0 <= w0,w1,w2 <= 1

        Vec3 intersection_point = ray.origin() + t * ray.direction();
        Vec3 local_point = intersection_point - m_c[0];

        float d00 = Vec3::dot(e1, e1);
        float d01 = Vec3::dot(e1, e2);
        float d11 = Vec3::dot(e2, e2);
        float d20 = Vec3::dot(local_point, e1);
        float d21 = Vec3::dot(local_point, e2);

        float denom = 1 / (d00 * d11 - d01 * d01);
        float v = (d11 * d20 - d01 * d21) * denom;
        float w = (d00 * d21 - d01 * d20) * denom;
        float u = 1.0f - v - w;

        // 4) Finally, if all conditions were met, then there is an intersection! :
        if (u >= 0 && v >= 0 && w >= 0 && u+v+w <= 1) {
            result.intersectionExists = true;
            result.t = t;
            result.intersection = intersection_point;
            result.normal = m_normal;
        }

        return result;
    }
};
#endif
