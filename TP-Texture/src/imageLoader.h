#ifndef IMAGELOADER_H
#define IMAGELOADER_H

#include <string>
#include <iostream>
#include <fstream>

// Source courtesy of J. Manson
// http://josiahmanson.com/prose/optimize_ppm/


namespace ppmLoader{
using namespace std;

void eat_comment(ifstream &f)
{
    char linebuf[1024];
    char ppp;
    while (ppp = f.peek(), ppp == '\n' || ppp == '\r')
        f.get();
    if (ppp == '#')
        f.getline(linebuf, 1023);
}


enum loadedFormat {
    rgb,
    rbg
};


void load_ppm( unsigned char * & pixels , unsigned int & w , unsigned int & h , const string &name , loadedFormat format = rgb)
{
    ifstream f(name.c_str(), ios::binary);
    if (f.fail())
    {
        cout << "Could not open file: " << name << endl;
        return;
    }

    // get type of file
    eat_comment(f);
    int mode = 0;
    string s;
    f >> s;
    if (s == "P3")
        mode = 3;
    else if (s == "P6")
        mode = 6;

    // get w
    eat_comment(f);
    f >> w;

    // get h
    eat_comment(f);
    f >> h;

    // get bits
    eat_comment(f);
    int bits = 0;
    f >> bits;

    // error checking
    if (mode != 3 && mode != 6)
    {
        cout << "Unsupported magic number" << endl;
        f.close();
        return;
    }
    if (w < 1)
    {
        cout << "Unsupported width: " << w << endl;
        f.close();
        return;
    }
    if (h < 1)
    {
        cout << "Unsupported height: " << h << endl;
        f.close();
        return;
    }
    if (bits < 1 || bits > 255)
    {
        cout << "Unsupported number of bits: " << bits << endl;
        f.close();
        return;
    }

    // load image data
    pixels = new unsigned char[3 * w * h];

    if (mode == 6)
    {
        f.get();
        f.read( (char*)&(pixels[0]) , 3 * w * h );
    }
    else if (mode == 3)
    {
        for (unsigned int i = 0; i < w * h; i++)
        {
            int v;
            f >> v;
            pixels[3*i + 0] = v;
            f >> v;
            pixels[3*i + 1] = v;
            f >> v;
            pixels[3*i + 2] = v;
        }
    }

    // close file
    f.close();

    if(format == rgb) {
        return;
    }

    if(format == rbg) {
        for (unsigned int i = 0; i < w * h; i++)
        {
            unsigned char r = pixels[3*i + 0] , g = pixels[3*i + 1] , b = pixels[3*i + 2];
            pixels[3*i + 0] = r;
            pixels[3*i + 1] = b;
            pixels[3*i + 2] = g;
        }
        return;
    }
}
}

#endif
