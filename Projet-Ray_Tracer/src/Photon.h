#ifndef PHOTON_H
#define PHOTON_H

#include "Vec3.h"

struct Photon {
    Vec3 position;      // Position of the photon
    Vec3 direction;     // Direction from which the photon came
    Vec3 power;         // Power (color) of the photon
    int bounceCount;    // Number of bounces from light source
    float flags;        // Bit flags for photon type (caustic, global, etc.)

    Photon(const Vec3& pos, const Vec3& dir, const Vec3& pow, int bounce = 0, float f = 0.0f)
        : position(pos), direction(dir), power(pow), 
          bounceCount(bounce), flags(f) {}

    Photon() 
        : position(Vec3(0,0,0)), direction(Vec3(0,0,0)), power(Vec3(0,0,0)), bounceCount(0), flags(0.0f) {}

    /*
    Photon(const Vec3& pos, const Vec3& dir, const Vec3& pow)
        : position(pos), direction(dir), power(pow) {}

    Photon() { position = Vec3(0, 0, 0); direction = Vec3(0, 0, 0); power = Vec3(0, 0, 0); }
    */

    // Helper methods for color operations
    void attenuate(const Vec3& surfaceColor) {
        power[0] *= surfaceColor[0];
        power[1] *= surfaceColor[1];
        power[2] *= surfaceColor[2];
    }

    void scaleEnergy(float factor) {
        power *= factor;
    }

    bool isSignificant() const {
        const float threshold = 0.001f;
        return (power[0] > threshold || 
                power[1] > threshold || 
                power[2] > threshold);
    }
    
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
