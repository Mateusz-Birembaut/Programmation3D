#ifndef SQUARE_H
#define SQUARE_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySquareIntersection{
    bool intersectionExists;
    float t;
    float u,v;
    Vec3 intersection;
    Vec3 normal;
};


class Square : public Mesh {
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;

    Square() : Mesh() {}
    Square(Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
           float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    void setQuad( Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
                  float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector , upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector*width;
        m_up_vector = m_up_vector*height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;                                      vertices[0].u = uMin; vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;                     vertices[1].u = uMax; vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;       vertices[2].u = uMax; vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;                        vertices[3].u = uMin; vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;


    }

    RaySquareIntersection intersect(const Ray &ray) const {
        RaySquareIntersection intersection;
        intersection.intersectionExists = false;

        Vec3 m_bottom_left = vertices[0].position;
        Vec3 m_right_vector = vertices[1].position - vertices[0].position;
        Vec3 m_up_vector = vertices[3].position - vertices[0].position;
        Vec3 m_normal = Vec3::cross(m_right_vector, m_up_vector);  
        m_normal.normalize();

        float D = Vec3::dot(m_bottom_left, m_normal);
        float denominator = Vec3::dot(ray.direction(), m_normal);

        if( denominator >= 0){ // si le determinant rayon.direction et plan.normal est nul, alors le rayon est parallele au plan donc pas d'intersection et si > 0, le rayon va "s'éloigner" donc pas d'intersection
            return intersection;
        }

        float t = (D - Vec3::dot(ray.origin(), m_normal))/ denominator;

        if( t <= 0){ // si l'intersection est derrière la caméra on s'en occupe pas
            return intersection;
        }

        Vec3 intersection_point = ray.origin() + t * ray.direction();
        // coordonnées du point d'intersection dans le repère local du carré
        Vec3 local_point = intersection_point - m_bottom_left;
        
        float x = Vec3::dot(local_point, m_right_vector) / (m_right_vector.length() *m_right_vector.length()); // donne le rapport de longueur entre le point projeté sur le vecteur droit et la longueur du vecteur droit 
        float y = Vec3::dot(local_point, m_up_vector) / (m_up_vector.length() *m_up_vector.length());
            
        if(x >= 0.0f && x <= 1 && y >= 0.0f && y <= 1){ // si on est dans le carré, on a une intersection
            intersection.intersectionExists = true;
            intersection.t = t;
            intersection.intersection = intersection_point;
            intersection.normal = m_normal;
            intersection.u = x;
            intersection.v = y;
        }

        return intersection;
    }

    
};
#endif // SQUARE_H
