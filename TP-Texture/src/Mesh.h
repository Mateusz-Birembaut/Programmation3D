#ifndef MESH_H
#define MESH_H


#include <vector>
#include <string>
#include "Vec3.h"

#include <GL/glut.h>


// -------------------------------------------
// Basic Mesh class
// -------------------------------------------

struct MeshVertex {
    inline MeshVertex () {}
    inline MeshVertex (const Vec3 & _p, const Vec3 & _n) : p (_p), n (_n) {}
    inline MeshVertex (const MeshVertex & v) : p (v.p), n (v.n) {}
    inline virtual ~MeshVertex () {}
    inline MeshVertex & operator = (const MeshVertex & v) {
        p = v.p;
        n = v.n;
        return (*this);
    }
    // membres :
    Vec3 p; // une position
    Vec3 n; // une normale
    Vec3 c; // une couleur

    float u;
    float v;

};

struct MeshTriangle {
    inline MeshTriangle () {
        corners[0] = corners[1] = corners[2] = 0;
    }
    inline MeshTriangle (const MeshTriangle & t) {
        corners[0] = t[0];   corners[1] = t[1];   corners[2] = t[2];
    }
    inline MeshTriangle (unsigned int v0, unsigned int v1, unsigned int v2) {
        corners[0] = v0;   corners[1] = v1;   corners[2] = v2;
    }
    inline virtual ~MeshTriangle () {}
    inline MeshTriangle & operator = (const MeshTriangle & t) {
        corners[0] = t[0];   corners[1] = t[1];   corners[2] = t[2];
        return (*this);
    }

    unsigned int operator [] (unsigned int c) const { return corners[c]; }
    unsigned int & operator [] (unsigned int c) { return corners[c]; }

private:
    // membres :
    unsigned int corners[3];
};




class Mesh {
public:
    std::vector<MeshVertex>   vertices;
    std::vector<MeshTriangle> triangles;
    std::vector<float>          positionArray;
    std::vector<unsigned int>   triangleArray;
    std::vector<float>          normalArray;
    std::vector<float>          colorArray;
    std::vector<float>          textureCoordArray;

    void loadOFF (const std::string & filename);
    void recomputeNormals ();
    void centerAndScaleToUnit ();
    void scaleUnit ();
    void buildVertexArray();
    void buildTriangleArray();
    void buildNormalArray();
    void buildColorArray();
    void setUnitSphere(int nX, int nY);
    void setSquare(int width, int height, float size, Vec3 center);

double square(double x) {return x*x;}


void setPseudoColors()
{
    for(unsigned int i=0; i<vertices.size();i++)
    {
        MeshVertex & v = vertices[i];
        double x = v.p[0];
        double y = v.p[1];
        double z = v.p[2];
        v.c[0] = ( square(square(y)) * square(sin(3* x )) );
        v.c[1] = ( square(cos(y)) );
        v.c[2] = ( square(y) * square(sin(3* z )) );
    }
}
void setUniformColor( Vec3 const & i_c )
{
    for(unsigned int i=0; i<vertices.size();i++)
    {
        MeshVertex & v = vertices[i];
        v.c = i_c;
    }
}

    void draw() const {
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glEnable(GL_COLOR_MATERIAL);
        glVertexPointer(3, GL_FLOAT, 3*sizeof(float), (GLvoid*)(&positionArray[0]));
        glNormalPointer(GL_FLOAT, 3*sizeof(float), (GLvoid*)(&normalArray[0]));
        glColorPointer(3, GL_FLOAT, 3*sizeof(float), (GLvoid*)(&colorArray[0]));
        glTexCoordPointer(2, GL_FLOAT, 2*sizeof(float), (GLvoid*)(&textureCoordArray[0]));
        
        glDrawElements(GL_TRIANGLES, triangleArray.size(), GL_UNSIGNED_INT, (GLvoid*)(&triangleArray[0]));
    }
};




#endif
