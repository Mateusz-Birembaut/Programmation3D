#ifndef GLOBALS_H
#define GLOBALS_H

//photon mapping
extern bool g_usePhotonMapping;
extern int g_photonCount;
extern float g_searchRadius;
extern int g_photonMaxBounces;

// camera
extern float g_camera_focalPlaneDistance;
extern float g_camera_apertureSize;

// ray tracing
extern int g_samplesPerPixel;
extern int g_shadowRays;
extern int g_rayMaxBounces;



#endif // GLOBALS_H