#version 120

varying vec4 P; // Per-vertex position
varying vec3 N; // Per-vertex normal
varying vec4 C; // Per-vertex color

uniform mat4 modelviewMatrix;
uniform mat4 projectionMatrix;
uniform mat3 normalMatrix;

// layout...
//layout(location = 1) in vec2 vertexUV;

//out sampler textureSampler;

void main(void) {
    P = gl_Vertex;
    N = gl_Normal;
    C = gl_Color;

    vec4 p = projectionMatrix * modelviewMatrix * P;
    gl_Position = p;
    //gl_MultiTexCoord0
    //textureSampler = ;


}
