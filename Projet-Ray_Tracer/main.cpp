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

#include "imgui/imgui.h"
#include "imgui/imgui_impl_glut.h"
#include "imgui/imgui_impl_opengl3.h"

#include <GL/glut.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <execution>
#include <chrono>
#include <algorithm>
#include <thread>
#include <vector>
#include <mutex>

#include <xmmintrin.h> // Inclure les intrinsics SIMD

#include "src/Scene.h"
#include "src/Vec3.h"
#include "src/Camera.h"
#include "src/matrixUtilities.h"
#include "src/imageLoader.h"
#include "src/Material.h"
#include "src/KdTree.h"
#include "src/Globals.h"

using namespace std;

// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------
static GLint window;
static unsigned int SCREENWIDTH = 480;
static unsigned int SCREENHEIGHT = 480;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX=0, lastY=0, lastZoom=0;
static unsigned int FPS = 0;
static bool fullScreen = false;

std::vector<Scene> scenes;
unsigned int selected_scene;

std::vector< std::pair< Vec3 , Vec3 > > rays;

void printUsage () {
    cerr << endl
         << "gMini: a minimal OpenGL/GLUT application" << endl
         << "for 3D graphics." << endl
         << "Author : Tamy Boubekeur (http://www.labri.fr/~boubek)" << endl << endl
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


// ------------------------------------

void initImGui() {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; 
    io.DisplaySize = ImVec2((float)SCREENWIDTH, (float)SCREENHEIGHT);
    ImGui::StyleColorsDark();
    ImGui_ImplGLUT_Init();
    ImGui_ImplOpenGL3_Init("#version 330");
}


void initLight () {
    GLfloat light_position[4] = {0.0, 1.5, 0.0, 1.0};
    GLfloat color[4] = { 1.0, 1.0, 1.0, 1.0};
    GLfloat ambient[4] = { 1.0, 1.0, 1.0, 1.0};

    glLightfv (GL_LIGHT1, GL_POSITION, light_position);
    glLightfv (GL_LIGHT1, GL_DIFFUSE, color);
    glLightfv (GL_LIGHT1, GL_SPECULAR, color);
    glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);
    glEnable (GL_LIGHT1);
    glEnable (GL_LIGHTING);
}

void init () {
    camera.resize (SCREENWIDTH, SCREENHEIGHT);
    initLight ();
    //glCullFace (GL_BACK);
    initImGui();
    glDisable (GL_CULL_FACE);
    glDepthFunc (GL_LESS);
    glEnable (GL_DEPTH_TEST);
    glClearColor (0.2f, 0.2f, 0.3f, 1.0f);
}


// ------------------------------------
// Replace the code of this 
// functions for cleaning memory, 
// closing sockets, etc.
// ------------------------------------

void clear () {

}

// ------------------------------------
// Replace the code of this 
// functions for alternative rendering.
// ------------------------------------

void addSpheresUI(Scene & scene) {
    std::vector<Sphere>& spheres = scene.getSpheres();
    if (spheres.size() == 0 ){
        return;
    }
    if (ImGui::CollapsingHeader("Spheres :")) {

        for (size_t i = 0; i < spheres.size(); i++){
            ImGui::PushID(static_cast<int>(i));
            if (ImGui::CollapsingHeader("Sphere ", i)) {
                Sphere & sphere = spheres[i];
                if (ImGui::SliderFloat("Radius", &sphere.m_radius, 0.0f, 10.0f))
                    scene.updateSphere(i);
                if (ImGui::SliderFloat3("Center", &sphere.m_center[0], -10.0f, 10.0f))
                    scene.updateSphere(i);
                if (sphere.material.type == Material_Glass){
                    ImGui::SliderFloat("Index of Refraction", &sphere.material.index_medium, 1.0f, 3.0f);
                }
                
                if (sphere.material.texture != nullptr){ 
                    ImGui::Text("Texture : %s", sphere.material.texture->name.c_str());
                    ImGui::InputFloat("Repeat x ", &sphere.material.t_uRepeat);
                    ImGui::InputFloat("Repeat y ", &sphere.material.t_vRepeat);
                    
                }
                if (sphere.material.normalMap != nullptr){ 
                    ImGui::Text("Normal Map : %s", sphere.material.normalMap->name.c_str());
                    ImGui::InputFloat("Repeat x ", &sphere.material.n_uRepeat);
                    ImGui::InputFloat("Repeat y ", &sphere.material.n_vRepeat);
                }
            }
            
            ImGui::PopID();
            ImGui::Separator();
        }
    }
}

void addSquaresUI(Scene & scene) {
    std::vector<Square>& squares = scene.getSquares();
    if (squares.size() == 0 ){
        return;
    }
    if (ImGui::CollapsingHeader("Squares :")) {
        for (size_t i = 0; i < squares.size(); i++){
            ImGui::PushID(static_cast<int>(i));
            if (ImGui::CollapsingHeader("Square ", i)) {
                Square & square = squares[i];
                if (square.material.type == Material_Glass){
                    ImGui::SliderFloat("Index of Refraction", &square.material.index_medium, 1.0f, 3.0f);
                }
                if (square.material.texture != nullptr){ 
                    if (ImGui::CollapsingHeader("Texture Settings")) {
                        ImGui::Text("Texture : %s", square.material.texture->name.c_str());
                        ImGui::InputFloat("Repeat x ", &square.material.t_uRepeat);
                        ImGui::InputFloat("Repeat y ", &square.material.t_vRepeat);
                    }
                }
                if (square.material.normalMap != nullptr){ 
                    if (ImGui::CollapsingHeader("Normal Map Settings")) {
                        ImGui::Text("Normal Map : %s", square.material.normalMap->name.c_str());
                        ImGui::InputFloat("Repeat x ", &square.material.n_uRepeat);
                        ImGui::InputFloat("Repeat y ", &square.material.n_vRepeat);
                    }
                }
            }
            
            ImGui::PopID();
            ImGui::Separator();
        }
    }
}

void addLightsUI(Scene & scene){
    std::vector<Light>& lights = scene.getLights();
    if (lights.size() == 0 ){
        return;
    }

    if (ImGui::CollapsingHeader("Lights :")) {
        for (size_t i = 0; i < lights.size(); i++){
            ImGui::PushID(static_cast<int>(i));
            if (ImGui::CollapsingHeader("Light ", i)) {
                Light & light = lights[i];
                ImGui::SliderFloat3("Position", &light.pos[0], -10.0f, 10.0f);
                ImGui::SliderFloat("Radius", &light.radius, 0.0f, 10.0f);
                ImGui::SliderFloat("Power Correction", &light.powerCorrection, 0.0f, 1000.0f);
                ImGui::SliderFloat3("Material", &light.material[0], 0.0f, 1.0f);
            }
            ImGui::PopID();
            ImGui::Separator();
        }
    }

}

void displayImGuiUI() {
    ImGui::Begin("Raytracer");

    ImGui::Text("Press 'r' to ray trace");
    ImGui::Separator(); 

    int temp_scene = static_cast<int>(selected_scene);
    ImGui::SliderInt("Scene", &temp_scene, 0, static_cast<int>(scenes.size() - 1));
    selected_scene = static_cast<unsigned int>(temp_scene);

    ImGui::InputInt("Samples per pixel", &g_samplesPerPixel);
    ImGui::InputInt("Shadow Rays per sample", &g_shadowRays);
    ImGui::InputInt("Max Ray Bounces", &g_rayMaxBounces);

    ImGui::Text("Objects dans la scene :");
    Scene & scene = scenes[selected_scene];

    addLightsUI(scene);
    addSpheresUI(scene);
    addSquaresUI(scene);


    ImGui::Separator();

    if (ImGui::CollapsingHeader("Camera Settings")) {
        ImGui::InputFloat("Focal Plane Distance", &g_camera_focalPlaneDistance);
        ImGui::InputFloat("Aperture Size", &g_camera_apertureSize);
    }

    ImGui::Separator();

    ImGui::Checkbox("Use Photon Mapping", &g_usePhotonMapping);
    if (g_usePhotonMapping) {
        if (ImGui::CollapsingHeader("Photon Mapping Settings")) {
            ImGui::InputInt("Photon Count", &g_photonCount, 0, 1000000);
            ImGui::SliderFloat("Search Radius", &g_searchRadius, 0.0f, 1.0f);
            ImGui::InputInt("Max Photon Bounces", &g_photonMaxBounces);
        }

    }


    ImGui::End();

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}


void draw () {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGLUT_NewFrame();
    ImGui::NewFrame();

    glEnable(GL_LIGHTING);

    scenes[selected_scene].draw();

    // draw rays : (for debug)
    //  std::cout << rays.size() << std::endl;
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    glLineWidth(6);
    glColor3f(1,0,0);
    glBegin(GL_LINES);
    for( unsigned int r = 0 ; r < rays.size() ; ++r ) {
        glVertex3f( rays[r].first[0],rays[r].first[1],rays[r].first[2] );
        glVertex3f( rays[r].second[0], rays[r].second[1], rays[r].second[2] );
    }

    glEnd();
    displayImGuiUI();

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
        sprintf (winTitle, "Raytracer - FPS: %d", FPS);
        glutSetWindowTitle (winTitle);
        lastTime = currentTime;
    }
    glutPostRedisplay ();
}


void ray_trace_block(int start_x, int end_x, int start_y, int end_y, float w, float h, unsigned int nsamples, std::vector<Vec3>& image, KdTreePhotonMap& kdTreePhotonMap, Vec3 pos) {
    float inv_w = 1.0f / w;
    float inv_h = 1.0f / h;
    float inv_nsamples = 1.0f / nsamples;
    for (int y = start_y; y < end_y; y++) {
        for (int x = start_x; x < end_x; x++) {
            Vec3 sum_color(0, 0, 0);
            for (unsigned int s = 0; s < nsamples; ++s) {
                float u = (x + dist05(rng)) * inv_w;
                float v = (y + dist05(rng)) * inv_h;
                Vec3 color = scenes[selected_scene].rayTrace(
                    kdTreePhotonMap,
                    depth_of_field_ray(u, v, camera.focalPlaneDistance, camera.apertureSize, pos)
                );
                sum_color += color;
            }
            image[x + y * w] = sum_color * inv_nsamples;
        }
    }
}


void ray_trace_from_camera() {
    int w = glutGet(GLUT_WINDOW_WIDTH)  ,   h = glutGet(GLUT_WINDOW_HEIGHT);
    std::cout << "Ray tracing a " << w << " x " << h << " image" << std::endl;
    camera.apply();
    Vec3 pos , dir;
    updateMatrices();
    pos = cameraSpaceToWorldSpace(Vec3(0,0,0));


    std::vector<Photon> photons;
    KdTreePhotonMap kdTreePhotonMap = KdTreePhotonMap(photons, 0);
    if (g_usePhotonMapping){
        scenes[selected_scene].photonMap(photons, g_photonCount); // x photons par source de lumière
        std::cout << "photons stocké  : " << photons.size() << std::endl;
        kdTreePhotonMap = KdTreePhotonMap(photons, 12);
    }
    photons.clear();

    camera.focalPlaneDistance = g_camera_focalPlaneDistance;
    camera.apertureSize = g_camera_apertureSize;


    int num_threads = std::thread::hardware_concurrency();
    int block_size_x = w / num_threads;
    int block_size_y = h / num_threads;

    unsigned int nsamples = g_samplesPerPixel;
    std::vector< Vec3 > image( w*h , Vec3(0,0,0) );
    auto start = std::chrono::high_resolution_clock::now();    

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        for (int j = 0; j < num_threads; ++j) {
            int start_x = i * block_size_x;
            int end_x = (i == num_threads - 1) ? w : start_x + block_size_x;
            int start_y = j * block_size_y;
            int end_y = (j == num_threads - 1) ? h : start_y + block_size_y;
            threads.emplace_back(ray_trace_block, start_x, end_x, start_y, end_y, w, h, nsamples, std::ref(image), std::ref(kdTreePhotonMap), pos);
        }
    }

    for (auto& thread : threads) {
        thread.join();
    } 

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Ray tracing completed in " << elapsed.count() << " seconds." << std::endl;
    std::cout << "\tDone" << std::endl;

    // ajouter post process ?

    std::string filename = "./rendu.ppm";
    ofstream f(filename.c_str(), ios::binary);
    if (f.fail()) {
        cout << "Could not open file: " << filename << endl;
        return;
    }
    f << "P3" << std::endl << w << " " << h << std::endl << 255 << std::endl;
    for (int i=0; i<w*h; i++)
        f << (int)(255.f*std::min<float>(1.f,image[i][0])) << " " << (int)(255.f*std::min<float>(1.f,image[i][1])) << " " << (int)(255.f*std::min<float>(1.f,image[i][2])) << " ";
    f << std::endl;
    f.close();

    char command[100];
    sprintf(command, "display rendu.ppm");
    if (system(command) == 0){
        std::cout << "Image displayed" << std::endl;
    }

}


void key (unsigned char keyPressed, int x, int y) {
    ImGui_ImplGLUT_KeyboardFunc(keyPressed, x, y);

    if (ImGui::GetIO().WantCaptureKeyboard) {
        return; 
    }

    Vec3 pos , dir;
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
        clear ();
        exit (0);
        break;
    case 'w':
        GLint polygonMode[2];
        glGetIntegerv(GL_POLYGON_MODE, polygonMode);
        if(polygonMode[0] != GL_FILL)
            glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
        else
            glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
        break;

    case 'r':
        camera.apply();
        rays.clear();
        ray_trace_from_camera();
        break;
    case '+':
        selected_scene++;
        if( selected_scene >= scenes.size() ) selected_scene = 0;
        break;
    default:
        printUsage ();
        break;
    }
    idle ();
}

void mouse (int button, int state, int x, int y) {
    ImGui_ImplGLUT_MouseFunc(button, state, x, y);
    if (ImGui::GetIO().WantCaptureMouse) {
        //std::cout << "ImGui wants to capture the mouse" << std::endl;
        return;
    }
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
    ImGui_ImplGLUT_MotionFunc(x, y); // Passer l'événement de mouvement de la souris à ImGui
    if (ImGui::GetIO().WantCaptureMouse) {
        return;
    }

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

void specialKey(int key, int x, int y) {
    ImGui_ImplGLUT_SpecialFunc(key, x, y);

    if (ImGui::GetIO().WantCaptureKeyboard) {
        return; 
    }

}

void specialKeyUp(int key, int x, int y) {
    ImGui_ImplGLUT_SpecialUpFunc(key, x, y);
}

void reshape(int w, int h) {
    camera.resize (w, h);
    SCREENWIDTH = w;
    SCREENHEIGHT = h;
    
    // Notifier ImGui de la nouvelle taille de la fenêtre
    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2((float)w, (float)h);
}

int main (int argc, char ** argv) {
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
    glutSpecialFunc(specialKey);
    glutReshapeFunc (reshape);
    glutMotionFunc (motion);
    glutMouseFunc (mouse);
    key ('?', 0, 0);

    camera.move(0., 0., -3.1);
    selected_scene=0;
    scenes.resize(6);
    scenes[0].setup_single_sphere();
    scenes[1].setup_single_square();
    scenes[2].setup_2_spheres();
    scenes[3].setup_2_planes();
    scenes[4].setup_cornell_box();
    scenes[5].setup_plan_2_spheres();

    glutMainLoop ();
    return EXIT_SUCCESS;
}

