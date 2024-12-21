#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include "Mesh.h"

#include <GL/glut.h>

class Scene {
    std::vector< Mesh > meshes;

public:
    Scene() {}

    void addMesh(std::string const & modelFilename) {
        meshes.resize( meshes.size() + 1 );
        Mesh & meshAjoute = meshes[ meshes.size() - 1 ];
        meshAjoute.loadOFF (modelFilename);
        meshAjoute.centerAndScaleToUnit ();
        meshAjoute.recomputeNormals ();
        meshAjoute.setUniformColor( Vec3(0.8,0.8,0.8) );

        meshAjoute.buildVertexArray();
        meshAjoute.buildTriangleArray();
        meshAjoute.buildNormalArray();
        meshAjoute.buildColorArray();
    }

    void addSphere( ) {
        meshes.resize( meshes.size() + 1 );
        Mesh & meshAjoute = meshes[ meshes.size() - 1 ];
        meshAjoute.setUnitSphere(20,20);
        meshAjoute.setPseudoColors();

        meshAjoute.buildVertexArray();
        meshAjoute.buildTriangleArray();
        meshAjoute.buildNormalArray();
        meshAjoute.buildColorArray();
    }


    void addSquare( int nX, int nY, float size, Vec3 center) {
        meshes.resize( meshes.size() + 1 );
        Mesh & meshAjoute = meshes[ meshes.size() - 1 ];
        meshAjoute.setSquare(nX+1,nY+1, size, center);
        meshAjoute.setPseudoColors();

        meshAjoute.buildVertexArray();
        meshAjoute.buildTriangleArray();
        meshAjoute.buildNormalArray();
        meshAjoute.buildColorArray();
    }

    void draw() const {
        // iterer sur l'ensemble des objets, et faire leur rendu.
        for( unsigned int mIt = 0 ; mIt < meshes.size() ; ++mIt ) {
            Mesh const & mesh = meshes[mIt];
            mesh.draw();
        }
    }
};



#endif
