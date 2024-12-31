#ifndef IMAGELOADER_H
#define IMAGELOADER_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "Vec3.h"
#include <algorithm>
#include "DataTypeEnum.h"

// Source courtesy of J. Manson
// http://josiahmanson.com/prose/optimize_ppm/


namespace ppmLoader{
using namespace std;
void eat_comment(ifstream &f);

struct RGB
{
    unsigned char r, g, b;
};



struct ImageRGB 
{
    int w, h;
    vector<RGB> data;
    string name;

    Vec3 samplePPMtoVec3(float u, float v , int repeatU, int repeatV, DATA_TYPE data_type) {
        u *= repeatU;
        v *= repeatV;

        u = u - std::floor(u);
        v = v - std::floor(v);

        if ( data_type == DATA_TYPE::NORMAL_MAP ) {
            v = 1.0f - v; // flip 
        }

        u = std::clamp(u, 0.0f, 1.0f);
        v = std::clamp(v, 0.0f, 1.0f);

        int x = static_cast<int>(u * (w - 1));
        int y = static_cast<int>(v * (h - 1));

        unsigned int index = y * w + x;

        if (index < 0 || index >= data.size()) {
            std::cerr << "Error: Texture coordinate out of bounds. u: " << u << ", v: " << v << ", index: " << index << std::endl;
            std::cerr << "Texture size: " << w << "x" << h << std::endl;
            return Vec3(0.0f, 0.0f, 0.0f);
        }

        RGB color = data[index];

        return Vec3(static_cast<float>(color.r) / 255.0f,
                    static_cast<float>(color.g) / 255.0f,
                    static_cast<float>(color.b) / 255.0f);
    }


};

void setName(ImageRGB & _img ,const string & _name);

void load_ppm(ImageRGB &img, const string &name);


enum loadedFormat {
    rgb,
    rbg
};


void load_ppm( unsigned char * & pixels , unsigned int & w , unsigned int & h , const string &name , loadedFormat format = rgb);


}


#endif
