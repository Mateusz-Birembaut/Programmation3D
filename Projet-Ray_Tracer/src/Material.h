#ifndef MATERIAL_H
#define MATERIAL_H

#include "imageLoader.h"
#include "Vec3.h"
#include <cmath>

#include <GL/glut.h>

enum MaterialType {
    Material_Diffuse_Blinn_Phong ,
    Material_Glass,
    Material_Mirror
};


struct Material {
    Vec3 ambient_material;
    Vec3 diffuse_material;
    Vec3 specular_material;
    double shininess;

    float index_medium;
    float transparency;

    MaterialType type;

    ppmLoader::ImageRGB * texture;
    float t_uRepeat;
    float t_vRepeat;

    ppmLoader::ImageRGB * normalMap;
    float n_uRepeat;
    float n_vRepeat;


    Material() {
        type = Material_Diffuse_Blinn_Phong;
        transparency = 0.0;
        index_medium = 1.0;
        ambient_material = Vec3(0., 0., 0.);
        texture = nullptr;
        normalMap = nullptr;
        shininess = 0.;
        n_uRepeat = 1.;
        n_vRepeat = 1.; 
        t_uRepeat = 1.;
        t_vRepeat = 1.; 
        
    }
};



#endif // MATERIAL_H
