#include "Mesh.h"
#include "Vec3.h"
#include <iostream>
#include <fstream>

void Mesh::loadOFF (const std::string & filename) {
    std::ifstream in (filename.c_str ());
    if (!in)
        exit (EXIT_FAILURE);
    std::string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
    vertices.resize (sizeV);
    triangles.resize (sizeT);
    for (unsigned int i = 0; i < sizeV; i++)
        in >> vertices[i].p;
    int s;
    for (unsigned int t = 0; t < sizeT; t++) {
        in >> s;
        for (unsigned int j = 0; j < 3; j++)
            in >> triangles[t][j];
    }
    in.close ();
}

void Mesh::recomputeNormals () {
    for (unsigned int i = 0; i < vertices.size (); i++)
        vertices[i].n = Vec3 (0.0, 0.0, 0.0);
    for (unsigned int t = 0; t < triangles.size (); t++) {
        Vec3 e01 = vertices[  triangles[t][1]  ].p -  vertices[  triangles[t][0]  ].p;
        Vec3 e02 = vertices[  triangles[t][2]  ].p -  vertices[  triangles[t][0]  ].p;
        Vec3 n = Vec3::cross (e01, e02);
        n.normalize ();
        for (unsigned int j = 0; j < 3; j++)
            vertices[  triangles[t][j]  ].n += n;
    }
    for (unsigned int i = 0; i < vertices.size (); i++)
        vertices[i].n.normalize ();
}

void Mesh::centerAndScaleToUnit () {
    Vec3 c(0,0,0);
    for  (unsigned int i = 0; i < vertices.size (); i++)
        c += vertices[i].p;
    c /= vertices.size ();
    float maxD = (vertices[0].p - c).length();
    for (unsigned int i = 0; i < vertices.size (); i++){
        float m = (vertices[i].p - c).length();
        if (m > maxD)
            maxD = m;
    }
    for  (unsigned int i = 0; i < vertices.size (); i++)
        vertices[i].p = (vertices[i].p - c) / maxD;
}

void Mesh::buildVertexArray()
{
    positionArray.clear();
    positionArray.reserve(3*vertices.size());
    for(unsigned int i=0; i<vertices.size();i++)
    {
        MeshVertex v = vertices[i];
        positionArray.push_back(v.p[0]);
        positionArray.push_back(v.p[1]);
        positionArray.push_back(v.p[2]);
    }
}

void Mesh::buildTriangleArray()
{
    triangleArray.clear();
    triangleArray.reserve(3*triangles.size());
    for(unsigned int i=0; i<triangles.size();i++)
    {
        MeshTriangle t = triangles[i];
        triangleArray.push_back(t[0]);
        triangleArray.push_back(t[1]);
        triangleArray.push_back(t[2]);
    }
}

void Mesh::buildNormalArray()
{
    normalArray.clear();
    normalArray.reserve(3*vertices.size());
    for(unsigned int i=0; i<vertices.size();i++)
    {
        MeshVertex v = vertices[i];
        normalArray.push_back(v.n[0]);
        normalArray.push_back(v.n[1]);
        normalArray.push_back(v.n[2]);
    }
}

void Mesh::buildColorArray()
{
    //Couleur aléatoire
    colorArray.clear();
    colorArray.reserve(3*vertices.size());
    for(unsigned int i=0; i<vertices.size();i++)
    {
        MeshVertex & v = vertices[i];
        colorArray.push_back( v.c[0] );
        colorArray.push_back( v.c[1] );
        colorArray.push_back( v.c[2] );
    }
}

void Mesh::setUnitSphere(int nX, int nY)
{
    vertices.clear();
    triangles.clear();
    float thetaStep = 2*M_PI/(nX-1);
    float phiStep = M_PI/(nY-1);
    for(int i=0; i<nX;i++)
    {
        for(int j=0;j<nY;j++)
        {
            float t = thetaStep*i;
            float p = phiStep*j - M_PI/2;

            Vec3 position(cos(t)*cos(p), sin(t)*cos(p), sin(p));
            Vec3 normal(position[0],position[1],position[2]);

            normal.normalize();
            vertices.push_back(MeshVertex(position, normal));
        }
    }
    for(int i=0; i<nX-1;i++)
    {
        for(int j=0;j<nY-1;j++)
        {
            triangles.push_back(MeshTriangle(i*nY+j, (i+1)*nY+j, (i+1)*nY+j+1));
            triangles.push_back(MeshTriangle(i*nY+j, (i+1)*nY+j+1, i*nY+j+1));
        }
    }

}


void Mesh::setSquare(int width, int height, float size, Vec3 center){
    vertices.clear();
    triangles.clear();

    //float stepWidth = float(1.0f/width-1);
    //float stepHeight = float(1.0f/height-1);

    float stepWidth = width;
    float stepHeight = height;

    for (int i = 0; i < width; i ++) {
        for (int j = 0; j < height; j ++) {
            
            Vec3 position( (i * stepWidth * size) + center[0], center[1], (j * stepHeight * size) + center[2]);
            Vec3 normal(0, 1 , 0);
            Vec3 color(1,1,1 );

            MeshVertex v = MeshVertex();

            v.c = color;
            v.n = normal;
            v.p = position;

            v.u = i * stepWidth;
            v.v = j * stepHeight;

            vertices.push_back(v);

        }
    }

    for(int i=0; i<width-1;i++){
        for(int j=0;j<height-1;j++){
            triangles.push_back(MeshTriangle((i+1)*height+j, i*height+j, (i+1)*height+j+1));
            triangles.push_back(MeshTriangle((i+1)*height+j+1, i*height+j , i*height+j+1));
        }
    }
}

