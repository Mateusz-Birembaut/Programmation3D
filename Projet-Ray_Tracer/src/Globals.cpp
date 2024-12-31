#include "Globals.h"

// Photon Mapping
bool g_usePhotonMapping = false;
int g_photonCount = 100000;
float g_searchRadius = 0.6f;
int g_photonMaxBounces = 5;

//Camera
float g_camera_focalPlaneDistance = 10.0f;
float g_camera_apertureSize = 0.0f;

//Ray tracing
int g_samplesPerPixel = 10;
int g_shadowRays = 5;
int g_rayMaxBounces = 3;
