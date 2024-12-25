#ifndef IMAGELOADER_H
#define IMAGELOADER_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "Vec3.h"
#include <algorithm>

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

    Vec3 sampleTextureAsVec3(float u, float v) {
        u = 1.0f - u; // flip u 
        std::clamp(u, 0.0f, 1.0f);
        std::clamp(v, 0.0f, 1.0f);

        int x = static_cast<int>(u * (w - 1));
        int y = static_cast<int>(v * (h - 1));

        unsigned int index = y * w + x;

        if (index < 0 || index >= data.size()) {
            std::cerr << "Error: Texture coordinate out of bounds. u: " << u << ", v: " << v << ", index: " << index << std::endl;
            return Vec3(0.0f, 0.0f, 0.0f);
        }

        RGB color = data[index];

        return Vec3(static_cast<float>(color.r) / 255.0f,
                    static_cast<float>(color.g) / 255.0f,
                    static_cast<float>(color.b) / 255.0f);
    }

};


void load_ppm(ImageRGB &img, const string &name);


enum loadedFormat {
    rgb,
    rbg
};


void load_ppm( unsigned char * & pixels , unsigned int & w , unsigned int & h , const string &name , loadedFormat format = rgb);


}


#endif
