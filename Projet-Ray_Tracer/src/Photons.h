#ifndef PHOTONS_H
#define PHOTONS_H

#include "Vec3.h"

struct Photon {
    Vec3 position;  // Position of the photon
    Vec3 direction; // Direction from which the photon came
    Vec3 power;     // Power (color) of the photon

    Photon(const Vec3& pos, const Vec3& dir, const Vec3& pow)
        : position(pos), direction(dir), power(pow) {}

    Photon() {}
};

#endif // PHOTONS_H