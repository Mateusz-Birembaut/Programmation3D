#version 330 core

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertices_position_modelspace;

//TODO create uniform transformations matrices Model View Projection
uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_projection;

// Values that stay constant for the whole mesh.

void main(){
        mat4 mvp = u_projection * u_view * u_model;
        // TODO : Output position of the vertex, in clip space : MVP * position
        gl_Position = mvp * vec4(vertices_position_modelspace,1)  ;

}

