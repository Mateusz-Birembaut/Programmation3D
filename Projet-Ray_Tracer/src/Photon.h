#ifndef PHOTON_H
#define PHOTON_H

#include "Vec3.h"

struct Photon {
    Vec3 position;  // Position of the photon
    Vec3 direction; // Direction from which the photon came
    Vec3 power;     // Power (color) of the photon

    Photon(const Vec3& pos, const Vec3& dir, const Vec3& pow)
        : position(pos), direction(dir), power(pow) {}

    Photon() { position = Vec3(0, 0, 0); direction = Vec3(0, 0, 0); power = Vec3(0, 0, 0); }
    
};

struct PhotonDistance {
    Photon photon;
    float distance;

    bool operator<(const PhotonDistance& other) const {
        return distance < other.distance;
    }

    PhotonDistance(const Photon& p, float d) : photon(p), distance(d) {}

};


#endif // PHOTONS_H
