#version 120

//const float attenuationConstant = 0.8;
//const float attenuationLinear = 0.8;
//const float attenuationQuadratic = 0.8;

uniform vec3 inputLightPosition;
uniform vec4 inputLightMaterial;
uniform bool lightIsInCamSpace;

uniform vec4 inputObjectAmbiantMaterial ;
uniform vec4 inputObjectDiffuseMaterial ;
uniform vec4 inputObjectSpecularMaterial ;
uniform float inputObjectShininess ;

//in sampler textureSampler;

uniform mat4 modelviewMatrix;
uniform mat4 projectionMatrix;
uniform mat3 normalMatrix;

varying vec4 P; // Interpolated fragment-wise position
varying vec3 N; // Interpolated fragment-wise normal
varying vec4 C; // Interpolated fragment-wise color

void main (void) {
    // PUT EVERY QUANTITY IN CAMERA SPACE!
    vec4 PModelview = modelviewMatrix * P;
    vec3 position = vec3 (PModelview.xyz / PModelview.w);
    vec3 normal = normalize(normalMatrix * N);
    vec3 directionToCamera = normalize(-position);
    vec3 lightPosition = inputLightPosition;
    if( !lightIsInCamSpace ) {
        vec4 LModelview = modelviewMatrix * vec4(lightPosition.xyz,1.0);
        lightPosition = vec3 (LModelview.xyz) / LModelview.w;
    }
    vec3 directionToLight = normalize(lightPosition - position);
    vec3 reflectedRayDirection = normalize( - directionToLight + 2.0*dot(directionToLight , normal) * normal );

    // ---------- Code to change -------------
    float distToLight = length(lightPosition - position);
    // float lightAttenuation = attenuationConstant + attenuationLinear/distToLight + attenuationQuadratic/(distToLight*distToLight);
    float lightAttenuation = 1.0;
    //vec4 color = C;
    vec4 color = C * lightAttenuation * ( (inputObjectAmbiantMaterial*inputLightMaterial)
                     + (inputObjectDiffuseMaterial*inputLightMaterial) * max(dot(directionToLight , normal),0.0)
                     + (inputObjectSpecularMaterial*inputLightMaterial) * pow( max(dot(reflectedRayDirection , directionToCamera) , 0.0) , inputObjectShininess )
                       );
    // ----------------------------------------
    color.w = 1.0;
    gl_FragColor = color;
}



