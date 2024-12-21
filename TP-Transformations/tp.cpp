// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>

// pour fmod (modulo avec flottants)
#include <cmath>

// Include GLEW
#include <GL/glew.h>

// Include GLM
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <GL/glut.h>

using namespace glm;
using Vec3 = glm::vec3;

using Mat4 = glm::mat4;

#include "src/shader.hpp"
#include "src/objloader.hpp"

// settings
const unsigned int SCR_WIDTH = 1333;
const unsigned int SCR_HEIGHT = 1000;

// camera
glm::vec3 camera_position   = glm::vec3(1.5f, 0.0f,  3.0f);
glm::vec3 camera_target = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 camera_up    = glm::vec3(0.0f, 1.0f,  0.0f);

// timing
float deltaTime = 0.1f;	// time between current frame and last frame
float lastFrame = 0.0f;

//rotation
float angle = 0.;
float zoom = 1.;

static GLint window;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX=0, lastY=0, lastZoom=0;
static bool fullScreen = false;

bool display_ex1 = true;
bool display_ex2 = false;
bool display_ex3 = false;

GLuint programID;
GLuint VertexArrayID;
GLuint vertexbuffer;
GLuint elementbuffer;
GLuint LightID;


std::vector<unsigned short> indices; //Triangles concaténés dans une liste
std::vector<std::vector<unsigned short> > triangles;
std::vector<glm::vec3> indexed_vertices;


glm::mat4 ViewMatrix;
glm::mat4 ProjectionMatrix;

glm::mat4 getViewMatrix(){
	return ViewMatrix;
}
glm::mat4 getProjectionMatrix(){
	return ProjectionMatrix;
}


// Initial position : on +Z
glm::vec3 position = glm::vec3( 0, 0, 0 );
// Initial horizontal angle : toward -Z
float horizontalAngle = 3.14f;
// Initial vertical angle : none
float verticalAngle = 0.0f;
// Initial Field of View
float initialFoV = 45.0f;

float speed = 3.0f; // 3 units / second
float mouseSpeed = 0.005f;


	// Right vector
glm::vec3 rightVector() {
    return glm::vec3(
		sin(horizontalAngle - 3.14f/2.0f),
		0,
		cos(horizontalAngle - 3.14f/2.0f)
	);
}

// Direction : Spherical coordinates to Cartesian coordinates conversion
glm::vec3 directionVector() {
    return glm::vec3(
        cos(verticalAngle) * sin(horizontalAngle),
        sin(verticalAngle),
        cos(verticalAngle) * cos(horizontalAngle)
    );
}

void computeMatricesFromInputs(float moveX, float moveY);
void initLight ();
void init ();
void draw ();
void display ();
void idle ();
void key (unsigned char keyPressed, int x, int y);
void mouse (int button, int state, int x, int y);
void motion (int x, int y);
void reshape(int w, int h);
int main (int argc, char ** argv);
void printMatrix(const glm::mat4& mat);

// ------------------------------------

void printMatrix(const glm::mat4& mat) {
    std::cout << mat[0][0] << " " << mat[1][0] << " " << mat[2][0] << " " << mat[3][0] << "\n" << mat[0][1] << " " << mat[1][1] << " " << mat[2][1] << " " << mat[3][1] << "\n" << mat[0][2] << " " << mat[1][2] << " " << mat[2][2] << " " << mat[3][2] << "\n" << mat[0][3] << " " << mat[1][3] << " " << mat[2][3] << " " << mat[3][3] << std::endl;
}

void initLight () {
    /*
    GLfloat light_position1[4] = {22.0f, 16.0f, 50.0f, 0.0f};
    GLfloat direction1[3] = {-52.0f,-16.0f,-50.0f};
    GLfloat color1[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat ambient[4] = {0.3f, 0.3f, 0.3f, 0.5f};
    */

    GLfloat light_position1[4] = {0.0f, 0.0f, 0.0f, 1.0f};
    GLfloat color1[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat ambient[4] = {0.3f, 0.3f, 0.3f, 0.5f};

    glLightfv (GL_LIGHT1, GL_POSITION, light_position1);
    //glLightfv (GL_LIGHT1, GL_SPOT_DIRECTION, direction1);
    glLightfv (GL_LIGHT1, GL_DIFFUSE, color1);
    glLightfv (GL_LIGHT1, GL_SPECULAR, color1);
    glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);
    glEnable (GL_LIGHT1);
    glEnable (GL_LIGHTING);
}

void init () {
    // camera.resize (SCREENWIDTH, SCREENHEIGHT);
    initLight ();
    //glCullFace (GL_BACK);
    glDisable (GL_CULL_FACE);
    glDepthFunc (GL_LESS);
    glEnable (GL_DEPTH_TEST);
    glClearColor (0.2f, 0.2f, 0.3f, 1.0f);
    glEnable(GL_COLOR_MATERIAL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // Initialize GLEW
    glewExperimental = true; // Needed for core profile
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        return;
    }

}


void DrawSphere(){
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(
                0,                  // attribute
                3,                  // size
                GL_FLOAT,           // type
                GL_FALSE,           // normalized?
                0,                  // stride
                (void*)0            // array buffer offset
                );

    // Index buffer
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);

    // Draw the triangles !
    glDrawElements(
                GL_TRIANGLES,      // mode
                indices.size(),    // count
                GL_UNSIGNED_SHORT,   // type
                (void*)0           // element array buffer offset
                );

    glDisableVertexAttribArray(0);
}

// ------------------------------------
// rendering.
// ------------------------------------


float multiplicateur_vitesse = 1.0f;

struct CorpsCeleste{
    float vitesseRotationOrbitale;    // vitesse de rotation autour du soleil
    float vitesseRotation;          // vitesse de rotation sur lui meme
    float taille;                // taille du corps 
    float angleRotation;         // angle de rotation sur lui meme
    float angleRotationOrbital;      // angle de rotation autour du soleil
    float inclinaisonAxiale;  // angle d'inclinaison de l'axe de rotation
    float distanceOrbite;   // distance par rapport au soleil
    float degresOrbite;     // degres de l'orbite (0 = orbite qui reste dans le plan xz)
    std::vector<CorpsCeleste> lunes; // liste des lunes

    CorpsCeleste(float taille, float vitesseRotationOrbitale, float vitesseRotation, float inclinaisonAxiale, float distanceOrbite, float degresOrbite) :
        vitesseRotationOrbitale(vitesseRotationOrbitale), 
        vitesseRotation(vitesseRotation),
        taille(taille), 
        inclinaisonAxiale(inclinaisonAxiale),
        distanceOrbite(distanceOrbite),
        degresOrbite(degresOrbite),
        angleRotation(0.0f),      
        angleRotationOrbital(0.0f)
    { }

    float getAngleRotation() {
        angleRotation = fmod(angleRotation + vitesseRotation * multiplicateur_vitesse, 360.0f);
        return angleRotation; 
    }

    float getAngleRotationOrbitale() { 
        angleRotationOrbital = fmod(angleRotationOrbital + vitesseRotationOrbitale * multiplicateur_vitesse, 360);
        return angleRotationOrbital; 
    }

    void addLune(const CorpsCeleste& lune) {
        lunes.push_back(lune);
    }

};


struct SystemeStellaire{
    std::vector<CorpsCeleste> planetes; // planetes du systeme
    CorpsCeleste etoile; // l'etoile du systeme

    SystemeStellaire(const CorpsCeleste& etoile) : 
        etoile(etoile), planetes() {
    }

    void addPlanete(const CorpsCeleste& planete) {
        planetes.push_back(planete);
    }

    void drawSystem() {
        glUseProgram(programID);

        // View matrix : camera/view transformation lookat() utiliser camera_position camera_target camera_up
        ViewMatrix = glm::lookAt(camera_position, camera_target, camera_up);

        // Projection matrix : 45 Field of View, 4:3 ratio, display range : 0.1 unit <-> 500 units
        ProjectionMatrix = glm::perspective(45.0f, 4.f / 3.f, 0.1f, 500.0f);

        //recup les loc des variables dans le shader
        GLint loc_transformations = glGetUniformLocation(programID, "u_model");
        GLint loc_ViewMatrix = glGetUniformLocation(programID, "u_view");
        GLint lov_ProjectionMatrix = glGetUniformLocation(programID, "u_projection");
        
        // on envoie les matrices view et projection au shader
        glUniformMatrix4fv(loc_ViewMatrix, 1, GL_FALSE, &ViewMatrix[0][0]);
        glUniformMatrix4fv(lov_ProjectionMatrix, 1, GL_FALSE, &ProjectionMatrix[0][0]);

        //matrice model pour l'étoile
        Mat4 model = Mat4(1.0);

        //rotation sur elle meme et angle si present
        model = glm::scale(model, glm::vec3(etoile.taille, etoile.taille, etoile.taille)); // changement taille
        model = glm::rotate(model, glm::radians(etoile.getAngleRotation()), glm::vec3(0, 1, 0)); // rotation sur y

        glUniformMatrix4fv(loc_transformations, 1, GL_FALSE, &model[0][0]);
        DrawSphere(); // dessine étoile

        // pour chaque planete
        for (CorpsCeleste& planete : planetes) {
            
            model = Mat4(1.0);//matrice model pour la planete actuelle

            model = glm::rotate(model, glm::radians(planete.degresOrbite), glm::vec3(0, 0, 1)); // degres entre plan xz et orbite planete
            model = glm::rotate(model, glm::radians(planete.getAngleRotationOrbitale()), glm::vec3(0, 1, 0)); // rotation autour du soleil
            model = glm::translate(model, glm::vec3(planete.distanceOrbite, 0, 0)); // translation selon distance avec le soleil

            for (CorpsCeleste& lune : planete.lunes) {
                Mat4 modelLune = model;

                modelLune = glm::rotate(modelLune, glm::radians(lune.degresOrbite), glm::vec3(0, 0, 1)); // degres entre orbit lune et orbite planete
                modelLune = glm::rotate(modelLune, glm::radians(lune.getAngleRotationOrbitale()), glm::vec3(0, 1, 0)); // orbite
                modelLune = glm::translate(modelLune, glm::vec3(  lune.distanceOrbite, 0, 0));  // déplacement selon distance avec la planete
                modelLune = glm::scale(modelLune, glm::vec3(lune.taille, lune.taille, lune.taille));  // changement taille
                modelLune = glm::rotate(modelLune, glm::radians(lune.getAngleRotation()), glm::vec3(0, 1 , 0)); // rotation sur elle meme

                glUniformMatrix4fv(loc_transformations, 1, GL_FALSE, &modelLune[0][0]); // envoie transfo pour placer la lune au shader
                DrawSphere(); // draw lune
            }

            //changement de la taille (apres la translation pour garder l'unité de distance cohérente entre toutes les planetes)
            model = glm::scale(model, glm::vec3(planete.taille, planete.taille, planete.taille)); 
            model = glm::rotate(model, glm::radians(planete.inclinaisonAxiale), glm::vec3(1, 0, 0));  // inclinaison axiale planete
            model = glm::rotate(model, glm::radians(planete.getAngleRotation()), glm::vec3(0, 1 , 0)); // rotation sur elle meme planete

            glUniformMatrix4fv(loc_transformations, 1, GL_FALSE, &model[0][0]);
            DrawSphere(); //planete
            
        }

    }
}; 


//CorpsCeleste(float taille, float vitesseRotationOrbitale, float vitesseRotation, float inclinaisonAxiale, float distanceOrbite, float degresOrbite) :

CorpsCeleste soleil = CorpsCeleste(5.f, 0.f, 0.000154f, 0.f, 0.0f, 0.0f); // Soleil, taille maximale

CorpsCeleste mercure = CorpsCeleste(1.f, 0.0000017f, 0.000017f, 0.03f, 4.f, 7.0f);    // Mercure
CorpsCeleste venus   = CorpsCeleste(1.8f, 0.0000011f, 0.000004f, 177.36f, 10.f, 3.39f);  // Vénus
CorpsCeleste terre   = CorpsCeleste(2.f, 0.0000114f, 0.0041667f, 23.44f, 16.f, 0.0f);    // Terre
CorpsCeleste lune    = CorpsCeleste(0.5f, 0.0000152f, 0.0000152f, 6.68f, 2.5f, 5.14f);   // Lune
CorpsCeleste mars    = CorpsCeleste(1.2f, 0.0000086f, 0.003521f, 25.19f, 20.f, 5.85f);   // Mars
CorpsCeleste jupiter = CorpsCeleste(4.f, 0.0000045f, 0.004545f, 3.13f, 26.f, 1.31f);  // Jupiter
CorpsCeleste saturne = CorpsCeleste(3.5f, 0.0000029f, 0.003684f, 26.73f, 32.f, 3.49f); // Saturne
CorpsCeleste uranus  = CorpsCeleste(2.5f, 0.0000015f, 0.001479f, 97.77f, 38.f, 0.77f); // Uranus
CorpsCeleste neptune = CorpsCeleste(2.4f, 0.0000011f, 0.001588f, 28.32f, 44.f, 5.77f);  // Neptune

// Mars
CorpsCeleste luneMarsPhobos = CorpsCeleste(0.22f*4, 0.0000011f, 0.0000011f, 9.4f, 5.0f, 0.3f); // Phobos
CorpsCeleste luneMarsDeimos = CorpsCeleste(0.12f*4, 0.000006f, 0.0000006f, 15.0f, 3.0f, 1.3f); // Deimos

// Jupiter
CorpsCeleste luneJupiterIo = CorpsCeleste(0.18f*4, 0.000092f, 0.0000892f, 10.0f, 1.8f*3, 1.8f); // Io
CorpsCeleste luneJupiterEurope = CorpsCeleste(0.16f*4, 0.0000489f, 0.0000489f, 18.5f, 3.5f*2, 3.5f); // Europe
CorpsCeleste luneJupiterGanymede = CorpsCeleste(0.26f*4, 0.0000148f, 0.0000148f, 7.0f, 0.9f*6, 0.9f); // Ganymède
CorpsCeleste luneJupiterCallisto = CorpsCeleste(0.24f*4, 0.0000158f, 0.0000158f, 16.7f, 1.0f*5, 1.0f); // Callisto

// Saturne
CorpsCeleste luneSaturneTitan = CorpsCeleste(0.23f*4, 0.0000134f, 0.0000134f, 10.5f, 5.0f, 5.0f); // Titan
CorpsCeleste luneSaturneRhea = CorpsCeleste(0.18f*4, 0.000041f, 0.000041f, 15.3f, 4.0f, 4.0f); // Rhea
CorpsCeleste luneSaturneIapetus = CorpsCeleste(0.14f*4, 0.000017f, 0.0000017f, 20.0f, 6.0f, 6.0f); // Iapetus
CorpsCeleste luneSaturneDione = CorpsCeleste(0.12f*4, 0.000025f, 0.0000025f, 25.0f, 7.0f, 7.0f); // Dione

// Uranus
CorpsCeleste luneUranusTitania = CorpsCeleste(0.18f*4, 0.000008f, 0.000008f, 20.0f, 8.0f, 8.0f); // Titania
CorpsCeleste luneUranusOberon = CorpsCeleste(0.16f*4, 0.000007f, 0.000007f, 25.0f, 9.0f, 9.0f); // Oberon
CorpsCeleste luneUranusAriel = CorpsCeleste(0.14f*4, 0.000005f, 0.000005f, 30.0f, 10.0f, 10.0f); // Ariel
CorpsCeleste luneUranusUmbriel = CorpsCeleste(0.13f*4, 0.000004f, 0.000004f, 35.0f, 11.0f, 11.0f); // Umbriel



// Neptune
CorpsCeleste luneNeptuneTriton = CorpsCeleste(0.14f*4, 0.000005f, 0.00005f, 18.0f, 4.5f, 4.5f); // Triton
CorpsCeleste luneNeptuneNereid = CorpsCeleste(0.2f*4, 0.0000001f, 0.0000001f, 18.5f, 6.0f, 6.0f); // Nereid



SystemeStellaire systeme = SystemeStellaire(soleil);


void drawEx2 () {
    glUseProgram(programID);
    // Model matrix : an identity matrix (model will be at the origin) then change
    Mat4 model = Mat4(1.0);

    // View matrix : camera/view transformation lookat() utiliser camera_position camera_target camera_up
    ViewMatrix = glm::lookAt(camera_position, camera_target, camera_up);

    // Projection matrix : 45 Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
    ProjectionMatrix = glm::perspective(45.0f, 4.f / 3.f, 0.1f, 100.0f);

    // Send our transformation to the currently bound shader,
    // in the "Model View Projection" to the shader uniforms
    GLint loc_transformations = glGetUniformLocation(programID, "u_model");
    GLint loc_ViewMatrix = glGetUniformLocation(programID, "u_view");
    GLint lov_ProjectionMatrix = glGetUniformLocation(programID, "u_projection");
    
    glUniformMatrix4fv(loc_ViewMatrix, 1,GL_FALSE, &ViewMatrix[0][0]);
    glUniformMatrix4fv(lov_ProjectionMatrix, 1,GL_FALSE, &ProjectionMatrix[0][0]);

    /*
   //ex 4.b
    // pour calculer l'angle : on fait u . v / norme(u) norme(v) = cos(theta)
    float produit_saclaire = glm::dot(glm::vec3(0., 1., 0.), glm::vec3(1., 1., 1.));
    float norme_u = glm::length(glm::vec3(0., 1., 0.));
    float norme_v = glm::length(glm::vec3(1., 1., 1.));

    float angle = acos(produit_saclaire / (norme_u * norme_v));

    glm::vec3 axe = glm::normalize(glm::cross(glm::vec3(0., 1., 0.), glm::vec3(1., 1., 1.))); 

    model =  glm::rotate(Mat4(1.0), angle, axe);
    */

    glUniformMatrix4fv(loc_transformations, 1,GL_FALSE, &model[0][0]);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(
                0,                  // attribute
                3,                  // size
                GL_FLOAT,           // type
                GL_FALSE,           // normalized?
                0,                  // stride
                (void*)0            // array buffer offset
                );

    // Index buffer
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);

    // Draw the triangles !
    glDrawElements(
                GL_TRIANGLES,      // mode
                indices.size(),    // count
                GL_UNSIGNED_SHORT,   // type
                (void*)0           // element array buffer offset
                );


    glDisableVertexAttribArray(0);
}

void drawChaise(){
    // 1rst attribute buffer : vertices
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(
                0,                  // attribute
                3,                  // size
                GL_FLOAT,           // type
                GL_FALSE,           // normalized?
                0,                  // stride
                (void*)0            // array buffer offset
                );

    // Index buffer
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);

    // Draw the triangles !
    glDrawElements(
                GL_TRIANGLES,      // mode
                indices.size(),    // count
                GL_UNSIGNED_SHORT,   // type
                (void*)0           // element array buffer offset
                );
}


float ex1_rotation = 0.0f;
void drawEx1(){
    glUseProgram(programID);
    // View matrix : camera/view transformation lookat() utiliser camera_position camera_target camera_up
    ViewMatrix = Mat4(1.0);

    // Projection matrix : 45 Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
    ProjectionMatrix = Mat4(1.0);

    // Send our transformation to the currently bound shader,
    // in the "Model View Projection" to the shader uniforms
    GLint loc_transformations = glGetUniformLocation(programID, "u_model");
    GLint loc_ViewMatrix = glGetUniformLocation(programID, "u_view");
    GLint lov_ProjectionMatrix = glGetUniformLocation(programID, "u_projection");
    
    glUniformMatrix4fv(loc_ViewMatrix, 1,GL_FALSE, &ViewMatrix[0][0]);
    glUniformMatrix4fv(lov_ProjectionMatrix, 1,GL_FALSE, &ProjectionMatrix[0][0]);

    // 1.A
    Mat4 transformations = Mat4(1.0);

    transformations[3][1] = -1.0; // sur le sol
    transformations[3][0] = -0.5; // un peu a gauche
    transformations[0][0] = 0.5; // scale x = 0.5
    transformations[1][1] = 0.5; // scale y = 0.5

    glUniformMatrix4fv(loc_transformations, 1,GL_FALSE, &transformations[0][0]);

    drawChaise();

    // Afficher une seconde chaise 1.B

    //ex1.b
    transformations = Mat4(1.0);

    transformations[3][1] = -1.0; // sur le sol
    transformations[3][0] = 0.5; // un peu a droite
    transformations[0][0] = -0.5; // scale x = -0.5 pour la retourner
    transformations[1][1] = 0.5; // scale y = 0.5

    glUniformMatrix4fv(loc_transformations, 1,GL_FALSE, &transformations[0][0]);

    drawChaise();

    // Afficher une troisieme chaise! 1.C

    //EX 1.3.C / D
    //matrice de translation pour avoir le centre de la chaise  en 0,0,0 
    Mat4 translation = Mat4(1.0);
    translation[3][0] = 0; 
    translation[3][1] = -0.5; 

    //matrice de rotation sur z (tourne sur elle meme au niveau du milieu de la chaise)
    Mat4 mat_rotation = Mat4(1.0);
    mat_rotation[0][0] = cos(glm::radians(ex1_rotation));
    mat_rotation[1][0] = -sin(glm::radians(ex1_rotation));
    mat_rotation[0][1] = sin(glm::radians(ex1_rotation));
    mat_rotation[1][1] = cos(glm::radians(ex1_rotation));

    //matrice de translation retour à la position initiale
    Mat4 translation_retour = Mat4(1.0);
    translation_retour[3][0] = 0; 
    translation_retour[3][1] = 0.5; 

    //translation_retour = glm::inverse(translation);

    // on met la chaise avec son centre de gravité au centre on rotationne et on remet la chaise à sa place
    transformations = translation_retour * mat_rotation * translation;

    //envoie les transformations au shader
    glUniformMatrix4fv(loc_transformations, 1,GL_FALSE, &transformations[0][0]);
    //dessine la chaise
    drawChaise();
    
    /*
    //ex 4.b
    // pour calculer l'angle : on utilise la formule du produit scalaire pour calculer l'angle entre deux vecteurs
    // u . v / norme(u) norme(v) = cos(theta)
    float produit_saclaire = glm::dot(glm::vec3(0., 1., 0.), glm::vec3(1., 1., 1.));
    float norme_u = glm::length(glm::vec3(0., 1., 0.));
    float norme_v = glm::length(glm::vec3(1., 1., 1.));

    float angle = acos(produit_saclaire / (norme_u * norme_v));

    glm::vec3 axe = glm::normalize(glm::cross(glm::vec3(0., 1., 0.), glm::vec3(1., 1., 1.))); 

    transformations =  glm::rotate(Mat4(1.0), angle, axe);

    glUniformMatrix4fv(loc_transformations, 1,GL_FALSE, &transformations[0][0]);
    */

    glDisableVertexAttribArray(0);
}

void draw () {
    if (display_ex1){
       drawEx1();
    }else if (display_ex2){
        drawEx2();
    }else if (display_ex3){ 
        systeme.drawSystem();
    }
}




void display () {
    glLoadIdentity ();
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // camera.apply ();
    draw ();
    glFlush ();
    glutSwapBuffers ();
}

void idle () {
    glutPostRedisplay ();
    float time = glutGet(GLUT_ELAPSED_TIME) / 1000.f;
    deltaTime = time - lastFrame;
    lastFrame = time;
}

void updateBuffers(){
        // mise a jour buffers
        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
        glBufferData(GL_ARRAY_BUFFER, indexed_vertices.size() * sizeof(glm::vec3), &indexed_vertices[0], GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned short), &indices[0] , GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void clearStoredMesh(){
    indexed_vertices.clear();
    indices.clear();
    triangles.clear();
}

void key (unsigned char keyPressed, int x, int y) {
    float cameraSpeed = 2.5 * deltaTime;
    std::string filename;
    glm::vec3 camera_direction = glm::normalize(camera_target - camera_position);
    glm::vec3 camera_right = glm::normalize(glm::cross(camera_direction, camera_up));
    switch (keyPressed) {
    case 'f':
        if (fullScreen == true) {
            glutReshapeWindow (SCR_WIDTH, SCR_HEIGHT);
            fullScreen = false;
        } else {
            glutFullScreen ();
            fullScreen = true;
        }
        break;

    case 's':
        camera_position -= cameraSpeed * camera_target ;
        break;

    case 'z':
        camera_position += cameraSpeed * camera_target ;
        break;
    case 'q':
        camera_position -= cameraSpeed * camera_right; 
        camera_target -= cameraSpeed * camera_right; 
        break;
    case 'd':
        camera_position += cameraSpeed * camera_right;
        camera_target += cameraSpeed * camera_right;  
        break;

    case 'a':
        ex1_rotation = fmod(ex1_rotation + 10, 360.0f);
        break;

    case 'e':
        ex1_rotation = fmod(ex1_rotation - 10, 360.0f);
        break;

    case '+':
        multiplicateur_vitesse *= 1.5f;
        break;

    case '-':
        multiplicateur_vitesse /= 1.5f;
        break;
    case '1':
        display_ex2 = false;
        display_ex3 = false;
        glDisable (GL_CULL_FACE);
        camera_position = glm::vec3(0.0f, 0.0f,  10.0f);
        clearStoredMesh();
        filename = "data/chair.off";
        loadOFF(filename, indexed_vertices, indices, triangles );
        updateBuffers();

        display_ex1 = !display_ex1;
        break;
    case '2':
        display_ex1 = false;
        display_ex3 = false;
        glEnable (GL_CULL_FACE);
        camera_position = glm::vec3(1.5f, 0.0f,  3.0f);
        filename = "data/suzanne.off";
        clearStoredMesh();
        loadOFF(filename, indexed_vertices, indices, triangles );
        updateBuffers();
        display_ex2 = !display_ex2;
        break;
    case '3':
        display_ex1 = false;
        display_ex2 = false;
        glEnable (GL_CULL_FACE);
        camera_position = glm::vec3(0.0f, 6.0f,  60.0f);
        filename = "data/sphere_2.off";
        clearStoredMesh();
        loadOFF(filename, indexed_vertices, indices, triangles );
        updateBuffers();
        display_ex3 = !display_ex3;
        break;

    default:
        break;
    }
    //TODO add translations
    idle ();
}

void specialKeys(int key, int x, int y) {
    if(key == GLUT_KEY_LEFT)
		position -= rightVector() * deltaTime * speed;
    else if(key == GLUT_KEY_RIGHT)
		position += rightVector() * deltaTime * speed;
    else if(key == GLUT_KEY_DOWN)
		position -= directionVector() * deltaTime * speed;
    else if(key == GLUT_KEY_UP)
        position += directionVector() * deltaTime * speed;
}

void mouse (int button, int state, int x, int y) {
    if (state == GLUT_UP) {
        mouseMovePressed = false;
        mouseRotatePressed = false;
        mouseZoomPressed = false;
    } else {
        if (button == GLUT_LEFT_BUTTON) {
            //camera.beginRotate (x, y);
            mouseMovePressed = false;
            mouseRotatePressed = true;
            mouseZoomPressed = false;
            lastX = x;
            lastY = y;
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
        computeMatricesFromInputs(x - lastX, y - lastY);
        lastX = x;
        lastY = y;
    }
    else if (mouseMovePressed == true) {
    }
    else if (mouseZoomPressed == true) {
    }
}

void computeMatricesFromInputs(float moveX, float moveY){
    std::cout << moveX << " " << moveY << std::endl;
	// Compute new orientation
	horizontalAngle += mouseSpeed * moveX / 10.f;
	verticalAngle   += mouseSpeed * moveY / 10.f;

	// Up vector
	glm::vec3 up = glm::cross( rightVector(), directionVector() );

	float FoV = initialFoV;

	// Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
	ProjectionMatrix = glm::perspective(glm::radians(FoV), 4.0f / 3.0f, 0.1f, 150.0f);

	// Camera matrix
	ViewMatrix       = glm::lookAt(
								camera_position,           // Camera is here
								camera_position + directionVector(), // and looks here : at the same position, plus "direction"
								up                  // Head is up (set to 0,-1,0 to look upside-down)
						   );
}


void reshape(int w, int h) {
    // camera.resize (w, h);
}

int main (int argc, char ** argv) {
    if (argc > 2) {
        exit (EXIT_FAILURE);
    }
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize (SCR_WIDTH, SCR_HEIGHT);
    window = glutCreateWindow ("TP HAI719I");

    init ();
    init ();
    glutIdleFunc(idle); // a
    glutDisplayFunc (display);
    glutKeyboardFunc (key);
    glutReshapeFunc (reshape);
    glutMotionFunc (motion);
    glutMouseFunc (mouse);
    glutSpecialFunc(specialKeys);
    key ('?', 0, 0);

    computeMatricesFromInputs(0.f, 0.f);

    //initialisation du systeme
    terre.addLune(lune);
    systeme.addPlanete(terre);

    mars.addLune(luneMarsPhobos);
    mars.addLune(luneMarsDeimos);
    systeme.addPlanete(mars);

    systeme.addPlanete(mercure);

    systeme.addPlanete(venus);

    jupiter.addLune(luneJupiterIo);
    jupiter.addLune(luneJupiterEurope);
    jupiter.addLune(luneJupiterGanymede);
    jupiter.addLune(luneJupiterCallisto);
    systeme.addPlanete(jupiter);

    saturne.addLune(luneSaturneTitan);
    saturne.addLune(luneSaturneRhea);
    saturne.addLune(luneSaturneIapetus);
    saturne.addLune(luneSaturneDione);
    systeme.addPlanete(saturne);

    uranus.addLune(luneUranusTitania);
    uranus.addLune(luneUranusOberon);
    uranus.addLune(luneUranusAriel);
    uranus.addLune(luneUranusUmbriel);
    systeme.addPlanete(uranus);

    neptune.addLune(luneNeptuneTriton);
    neptune.addLune(luneNeptuneNereid);
    systeme.addPlanete(neptune); 

    glGenVertexArrays(1, &VertexArrayID);
    glBindVertexArray(VertexArrayID);

    // Create and compile our GLSL program from the shaders
    programID = LoadShaders( "vertex_shader.glsl", "fragment_shader.glsl" );

    //Chargement du fichier de maillage
    std::string filename("data/chair.off");
    //std::string filename("data/sphere_2.off");
    loadOFF(filename, indexed_vertices, indices, triangles );

    // Load it into a VBO

    glGenBuffers(1, &vertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, indexed_vertices.size() * sizeof(glm::vec3), &indexed_vertices[0], GL_STATIC_DRAW);

    // Generate a buffer for the indices as well
    glGenBuffers(1, &elementbuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned short), &indices[0] , GL_STATIC_DRAW);

    // Get a handle for our "LightPosition" uniform
    glUseProgram(programID);
    LightID = glGetUniformLocation(programID, "LightPosition_worldspace");

    glutMainLoop ();

    // Cleanup VBO and shader
    glDeleteBuffers(1, &vertexbuffer);
    glDeleteBuffers(1, &elementbuffer);
    glDeleteProgram(programID);
    glDeleteVertexArrays(1, &VertexArrayID);


    return EXIT_SUCCESS;
}
