// -------------------------------------------
// gMini : a minimal OpenGL/GLUT application
// for 3D graphics.
// Copyright (C) 2006-2008 Tamy Boubekeur
// All rights reserved.
// -------------------------------------------

// -------------------------------------------
// Disclaimer: this code is dirty in the
// meaning that there is no attention paid to
// proper class attribute access, memory
// management or optimisation of any kind. It
// is designed for quick-and-dirty testing
// purpose.
// -------------------------------------------

#include "src/GLProgram.h"
#include "src/Exception.h"

#include <GL/glext.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include "src/Mesh.h"
#include "src/Vec3.h"
#include "src/Camera.h"
#include "src/Scene.h"
#include <GL/glut.h>

using namespace std;

#include "src/imageLoader.h"


// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------

static GLint window;
static unsigned int SCREENWIDTH = 640;
static unsigned int SCREENHEIGHT = 480;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX=0, lastY=0, lastZoom=0;
static unsigned int FPS = 0;
static bool fullScreen = false;

Scene scene;
GLProgram * glProgram;
GLProgram * glProgramTextures;


struct Texture{
    unsigned char * values;
    unsigned int width;
    unsigned int height;
    GLuint handleIndex;
    GLuint bindingIndex;
    GLuint format = GL_RGBA; // could be GL_RGB
};


Texture checkerBoardColorTexture;
Texture normalMapTexture;


void printUsage () {
    cerr << endl
         << "Usage : ./gmini [<file.off>]" << endl
         << "Keyboard commands" << endl
         << "------------------" << endl
         << " ?: Print help" << endl
         << " w: Toggle Wireframe Mode" << endl
         << " g: Toggle Gouraud Shading Mode" << endl
         << " f: Toggle full screen mode" << endl
         << " <drag>+<left button>: rotate model" << endl
         << " <drag>+<right button>: move model" << endl
         << " <drag>+<middle button>: zoom" << endl
         << " q, <esc>: Quit" << endl << endl;
}

void usage () {
    printUsage ();
    exit (EXIT_FAILURE);
}

void createCheckerBoardImage( ) {
    checkerBoardColorTexture.width = checkerBoardColorTexture.height = 256; // for example
    checkerBoardColorTexture.values = new unsigned char[ 4 * checkerBoardColorTexture.width * checkerBoardColorTexture.height ];
    for(unsigned int i = 0 ; i < checkerBoardColorTexture.width ; ++i) {
        for(unsigned int j=0 ; j < checkerBoardColorTexture.height ; ++j) {
            char value = 0;
            if(((2*i)/checkerBoardColorTexture.width + (2*j)/checkerBoardColorTexture.height) %2) {
                value = 255u;
            }
            checkerBoardColorTexture.values[4 * (i*checkerBoardColorTexture.height+j) +0] = value; // red
            checkerBoardColorTexture.values[4 * (i*checkerBoardColorTexture.height+j) +1] = value; // green
            checkerBoardColorTexture.values[4 * (i*checkerBoardColorTexture.height+j) +2] = value; // blue
            checkerBoardColorTexture.values[4 * (i*checkerBoardColorTexture.height+j) +3] = 255;   // alpha ( = 1 )
        }
    }
}

void init () {
    glewExperimental = GL_TRUE;
    glewInit (); // init glew, which takes in charges the modern OpenGL calls (v>1.2, shaders, etc)
    camera.resize (SCREENWIDTH, SCREENHEIGHT);
    glCullFace (GL_BACK);
    glEnable (GL_CULL_FACE);
    glDepthFunc (GL_LESS);
    glEnable (GL_DEPTH_TEST);
    glClearColor (0.2f, 0.2f, 0.3f, 1.0f);

    // load images:
    createCheckerBoardImage();
    ppmLoader::load_ppm(normalMapTexture.values , normalMapTexture.width , normalMapTexture.height , "./img/normalMaps/n3.ppm" , ppmLoader::rbg);
    normalMapTexture.format = GL_RGB; // only three channels in the normal map



    //----------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------------------------//
    //-------------------------------------         GENERATE YOUR TEXTURES HERE !       ------------------------------------------//
    glEnable(GL_TEXTURE_2D);


    GLuint Texture;
    glGenTextures(1, &Texture);

    glBindTexture(GL_TEXTURE_2D, Texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, normalMapTexture.width, normalMapTexture.height, 0, GL_RGB, GL_UNSIGNED_BYTE, normalMapTexture.values);
    glBindTexture(GL_TEXTURE_2D, 0);


    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, Texture);
    glUniform1i(Texture, 0);

    glDisable(GL_TEXTURE_2D);
    //----------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------------------------------------------------//


    //scene.addSquare(1,1);



    // Load first shader
    try {
        glProgram = GLProgram::genVFProgram ("Simple GL Program", "./src/shader.vert", "./src/shader.frag");
        // Load and compile pair of shaders
        glProgram->use ();



        glProgram->stop ();

        // Activate the shader program
    } catch (Exception & e) {
        cerr << e.msg () << endl;
    }
}


void draw () {
    GLfloat mvmatrix[16];
    GLfloat nmatrix[9];
    GLfloat pjmatrix[16];

    glMatrixMode(GL_PROJECTION);
    glGetFloatv (GL_PROJECTION_MATRIX,pjmatrix);
    glMatrixMode(GL_MODELVIEW);
    glGetFloatv (GL_MODELVIEW_MATRIX, mvmatrix);
    for( int i = 0 ; i < 3 ; ++i ) {
        for( int j = 0 ; j < 3 ; ++j ) {
            nmatrix[3*i+j] = mvmatrix[4*i+j];
        }
    }


    // set openGL parameters:
    glProgram->use();
    {
        glProgram->setUniform3f("inputLightPosition" , 1.0, 1.0, 1.0);
        glProgram->setUniform4f("inputLightMaterial" , 1.0, 1.0, 1.0 , 1.0);
        glProgram->setUniform1i("lightIsInCamSpace" , 0);

        glProgram->setUniform4f("inputObjectAmbiantMaterial" , 0.1, 0.2, 0.1 , 1.0);
        glProgram->setUniform4f("inputObjectDiffuseMaterial" , 0.3, 0.3, 0.4 , 1.0);
        glProgram->setUniform4f("inputObjectSpecularMaterial" , 1.0, 1.0, 1.0 , 1.0);
        glProgram->setUniform1f("inputObjectShininess" , 50.0);

        glProgram->setUniformMatrix4fv("modelviewMatrix",mvmatrix);
        glProgram->setUniformMatrix4fv("projectionMatrix",pjmatrix);
        glProgram->setUniformMatrix3fv("normalMatrix",nmatrix);

        //glProgram->setUniform


    }
    //{
        //glBegin(GL_QUADS);
        //glTexCoord2f(-10,10); glVertex3f(-10, -1.5, 10);
        //glTexCoord2f(10,10); glVertex3f(10, -1.5, 10);
        //glTexCoord2f(10,-10); glVertex3f(10, -1.5, -10);
        //glTexCoord2f(-10,-10); glVertex3f(-10, -1.5, -10);

        //glEnd();
    //}

    Vec3 center = Vec3(-5,-1,-5);
    scene.addSquare(1, 1, 10.,center);
    scene.draw();
    glProgram->stop();


}

void display () {
    glLoadIdentity ();
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply ();
    draw ();
    glFlush ();
    glutSwapBuffers ();
}

void idle () {
    static float lastTime = glutGet ((GLenum)GLUT_ELAPSED_TIME);
    static unsigned int counter = 0;
    counter++;
    float currentTime = glutGet ((GLenum)GLUT_ELAPSED_TIME);
    if (currentTime - lastTime >= 1000.0f) {
        FPS = counter;
        counter = 0;
        static char winTitle [64];
        sprintf (winTitle, "gmini - FPS: %d", FPS);
        glutSetWindowTitle (winTitle);
        lastTime = currentTime;
    }
    glutPostRedisplay ();
}

void key (unsigned char keyPressed, int x, int y) {
    switch (keyPressed) {
    case 'f':
        if (fullScreen == true) {
            glutReshapeWindow (SCREENWIDTH, SCREENHEIGHT);
            fullScreen = false;
        } else {
            glutFullScreen ();
            fullScreen = true;
        }
        break;
    case 'q':
    case 27:
        exit (0);
        break;
    case 'w':
        glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
        break;
    case 'g':
        glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
        break;

    default:
        printUsage ();
        break;
    }
    idle ();
}

void mouse (int button, int state, int x, int y) {
    if (state == GLUT_UP) {
        mouseMovePressed = false;
        mouseRotatePressed = false;
        mouseZoomPressed = false;
    } else {
        if (button == GLUT_LEFT_BUTTON) {
            camera.beginRotate (x, y);
            mouseMovePressed = false;
            mouseRotatePressed = true;
            mouseZoomPressed = false;
        } else if (button == GLUT_RIGHT_BUTTON) {
            lastX = x;
            lastY = y;
            mouseMovePressed = true;
            mouseRotatePressed = false;
            mouseZoomPressed = false;
        } else if (button == GLUT_MIDDLE_BUTTON) {
            if (mouseZoomPressed == false) {
                lastZoom = y;
                mouseMovePressed = false;
                mouseRotatePressed = false;
                mouseZoomPressed = true;
            }
        }
    }
    idle ();
}

void motion (int x, int y) {
    if (mouseRotatePressed == true) {
        camera.rotate (x, y);
    }
    else if (mouseMovePressed == true) {
        camera.move ((x-lastX)/static_cast<float>(SCREENWIDTH), (lastY-y)/static_cast<float>(SCREENHEIGHT), 0.0);
        lastX = x;
        lastY = y;
    }
    else if (mouseZoomPressed == true) {
        camera.zoom (float (y-lastZoom)/SCREENHEIGHT);
        lastZoom = y;
    }
}


void reshape(int w, int h) {
    camera.resize (w, h);
}


int main (int argc, char ** argv) {
    std::cout << "start" << std::endl;
    if (argc > 2) {
        printUsage ();
        exit (EXIT_FAILURE);
    }
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize (SCREENWIDTH, SCREENHEIGHT);
    window = glutCreateWindow ("gMini");

    init ();
    glutIdleFunc (idle);
    glutDisplayFunc (display);
    glutKeyboardFunc (key);
    glutReshapeFunc (reshape);
    glutMotionFunc (motion);
    glutMouseFunc (mouse);
    key ('?', 0, 0);

    scene.addMesh (argc == 2 ? argv[1] : "models/monkey.off");

    glutMainLoop ();
    return EXIT_SUCCESS;
}

