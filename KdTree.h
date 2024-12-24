#ifndef KDTREE_H
#define KDTREE_H


#include <vector>
#include "BoundingBox.h"
#include "Sphere.h"
#include "Square.h"




class KdTreeNode {
public:

    KdTreeNode* left;
    KdTreeNode* right;

    int dimSplit; // x y z = 0 1 2
    float splitDistance;

    bool isLeaf = false;
    std::vector<Triangle> triangles; // modif pour avoir des adresses de triangles
    BoundingBox box;

    KdTreeNode() : left(nullptr), right(nullptr) {}

    KdTreeNode* KdTreeBuild(const std::vector<Triangle>& _triangles, int _depth, BoundingBox parentBox = BoundingBox()){
        if (_triangles.size() == 0){
            return nullptr;
        }

        KdTreeNode* node = new KdTreeNode();
        
        if (parentBox.min == Vec3(0,0,0) && parentBox.max == Vec3(0,0,0)) {
            node->box = BoundingBox().meshBoundingBox(_triangles);
        } else {
            node->box = parentBox;
        }
        //std::cout << "Node box: " << node->box.min << " " << node->box.max << std::endl;

        if (_depth == 0 || _triangles.size() <= 2) {
            node->isLeaf = true;
            for (int i = 0; i < _triangles.size(); i++){
                node->triangles.push_back(_triangles[i]);
            }
            return node;
        }

        Vec3 diff = node->box.max - node->box.min;
        node->dimSplit = diff.maxDimension();
        node->splitDistance = node->box.min[node->dimSplit] + diff[node->dimSplit] / 2;

        BoundingBox leftBox = node->box;
        BoundingBox rightBox = node->box;
        leftBox.max[node->dimSplit] = node->splitDistance;
        rightBox.min[node->dimSplit] = node->splitDistance;
        //std::cout << "Left box: " << leftBox.min << " " << leftBox.max << std::endl;
        //std::cout << "Right box: " << rightBox.min << " " << rightBox.max << std::endl;

        std::vector<Triangle> leftTriangles;
        std::vector<Triangle> rightTriangles;

        for (const Triangle& triangle : _triangles){ 
            BoundingBox triangleBox = BoundingBox().triangleBoundingBox(triangle);
            
            if (triangleBox.max[node->dimSplit] < node->splitDistance){
                leftTriangles.push_back(triangle);
            } else if (triangleBox.min[node->dimSplit] > node->splitDistance){
                rightTriangles.push_back(triangle);
            } else {
                leftTriangles.push_back(triangle);
                rightTriangles.push_back(triangle);
            }
        }

        //std::cout << "Left triangles: " << leftTriangles.size() << " Right triangles: " << rightTriangles.size() << std::endl;

        node->left = KdTreeBuild(leftTriangles, _depth - 1, leftBox);
        node->right = KdTreeBuild(rightTriangles, _depth - 1, rightBox);

        return node;
    }

    void traverse(Ray const & ray, RayTriangleIntersection & intersection) const { 
        float tEntree = box.intersect(ray);
        if (tEntree != INFINITY) {
            if (left == nullptr && right == nullptr) {
                //std::cout << "Leaf node with " << triangles.size() << " triangles." << std::endl;
                for (const Triangle& triangle : triangles) {
                    //std::cout << "Triangle vertices: " << triangle.getVertices()[0] << " " << triangle.getVertices()[1]<< " " << triangle.getVertices()[2] << std::endl;
                    RayTriangleIntersection tempIntersection = triangle.getIntersection(ray);
                    if (tempIntersection.intersectionExists && tempIntersection.t < intersection.t && tempIntersection.t > 0.001) {
                        //std::cout << "Intersection found at distance " << tempIntersection.t << std::endl;
                        //std::cout << " Distance avec BB: " << tEntree << std::endl;
                        intersection = tempIntersection;
                    }
                }
            } else {
                if (left != nullptr ){
                    left->traverse(ray, intersection);
                }
                if (right != nullptr){
                    right->traverse(ray, intersection);
                }
            }
        }
        
    }


    void testKdTreeBuild() {
        std::vector<Triangle> triangles;
        triangles.push_back(Triangle(Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0)));
        triangles.push_back(Triangle(Vec3(1, 1, 1), Vec3(2, 1, 1), Vec3(1, 2, 1)));
        triangles.push_back(Triangle(Vec3(2, 2, 2), Vec3(3, 2, 2), Vec3(2, 3, 2)));
        //triangles.push_back(Triangle(Vec3(12, 12, 12), Vec3(11, 11, 11), Vec3(10, 10, 10)));


        KdTreeNode* root = nullptr;
        root = root->KdTreeBuild(triangles, 3);

        if (root->isLeaf) {
            std::cout << "Root is a leaf node with " << root->triangles.size() << " triangles." << std::endl;
        } else {
            std::cout << "Root is not a leaf node." << std::endl;
        }
    }

    void testTraverse() {
        std::vector<Triangle> testTriangles;
        testTriangles.push_back(Triangle(Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0)));
        testTriangles.push_back(Triangle(Vec3(1, 1, 1), Vec3(2, 1, 1), Vec3(1, 2, 1)));

        KdTreeNode* root = nullptr;
        root = root->KdTreeBuild(testTriangles, 3);

        Ray ray(Vec3(0, 0, -1), Vec3(0, 0, 1));
        RayTriangleIntersection intersection;

        root->traverse(ray, intersection);
        if (intersection.intersectionExists)
        {
            std::cout << "Intersection found at distance " << intersection.t << std::endl;
        }else{
            std::cout << "No intersection found." << std::endl;
        }
    
    }   


    void testKdTreeBuildComplex() {
        std::vector<Triangle> triangles;

        triangles.push_back(Triangle(Vec3(-1, -1, 1), Vec3(-0.5, -1, 0), Vec3(-1, -0.5, 0)));
        triangles.push_back(Triangle(Vec3(0, 0, 0),   Vec3(1, 0, 0),    Vec3(0, 1, 0)));
        triangles.push_back(Triangle(Vec3(2, 2, 0),   Vec3(3, 2, 0),    Vec3(2, 3, 0)));
        triangles.push_back(Triangle(Vec3(5, 5, 5),   Vec3(4, 6, 5),    Vec3(5, 6, 5)));
        triangles.push_back(Triangle(Vec3(10, 10, -2),Vec3(12, 10, -2), Vec3(10, 12, -2)));
        triangles.push_back(Triangle(Vec3(-2, 2, 2),  Vec3(-3, 2, 2),   Vec3(-2, 3, 2)));

        KdTreeNode* root = nullptr;

        root = root->KdTreeBuild(triangles, 4);

        if (!root) {
            std::cout << "Error: root is null." << std::endl;
            return;
        }

        if (root->isLeaf) {
            std::cout << "Root is a leaf node with " << root->triangles.size() << " triangles." << std::endl;
        } else {
            std::cout << "Root is not a leaf node." << std::endl;
        }
        
        // Testez un rayon traversant plusieurs zones
        Ray ray(Vec3(-10, -10, 0), Vec3(1, 1, 0));
        RayTriangleIntersection intersection;
        intersection.intersectionExists = false;
        intersection.t = std::numeric_limits<float>::max();

        root->traverse(ray, intersection);
        if (intersection.intersectionExists)
        {
            std::cout << "Intersection found: " << intersection.intersectionExists 
                << " at distance " << intersection.t << std::endl;
        }else {
            std::cout << "No intersection found." << std::endl;
        }
    } 

    void drawBoundingBoxesHelper(const KdTreeNode* node) const {
        if (!node) return;

        // Draw the triangles if the node is a leaf
        if (node->left == nullptr && node->right == nullptr) {
            // Set color for better visibility
            glColor3f(1.0f, 1.0f, 1.0f);  // Red color
            
            // Enable line width
            glLineWidth(1.0f);

            glBegin(GL_LINES);


            // Draw the edges of the bounding box
            glVertex3f(node->box.min[0], node->box.min[1], node->box.min[2]);
            glVertex3f(node->box.max[0], node->box.min[1], node->box.min[2]);

            glVertex3f(node->box.min[0], node->box.min[1], node->box.min[2]);
            glVertex3f(node->box.min[0], node->box.max[1], node->box.min[2]);

            glVertex3f(node->box.min[0], node->box.min[1], node->box.min[2]);
            glVertex3f(node->box.min[0], node->box.min[1], node->box.max[2]);

            glVertex3f(node->box.max[0], node->box.max[1], node->box.max[2]);
            glVertex3f(node->box.min[0], node->box.max[1], node->box.max[2]);

            glVertex3f(node->box.max[0], node->box.max[1], node->box.max[2]);
            glVertex3f(node->box.max[0], node->box.min[1], node->box.max[2]);

            glVertex3f(node->box.max[0], node->box.max[1], node->box.max[2]);
            glVertex3f(node->box.max[0], node->box.max[1], node->box.min[2]);

            glVertex3f(node->box.min[0], node->box.max[1], node->box.min[2]);
            glVertex3f(node->box.max[0], node->box.max[1], node->box.min[2]);

            glVertex3f(node->box.min[0], node->box.max[1], node->box.min[2]);
            glVertex3f(node->box.min[0], node->box.max[1], node->box.max[2]);

            glVertex3f(node->box.max[0], node->box.min[1], node->box.min[2]);
            glVertex3f(node->box.max[0], node->box.max[1], node->box.min[2]);

            glVertex3f(node->box.max[0], node->box.min[1], node->box.min[2]);
            glVertex3f(node->box.max[0], node->box.min[1], node->box.max[2]);

            glVertex3f(node->box.min[0], node->box.min[1], node->box.max[2]);
            glVertex3f(node->box.max[0], node->box.min[1], node->box.max[2]);

            glVertex3f(node->box.min[0], node->box.min[1], node->box.max[2]);
            glVertex3f(node->box.min[0], node->box.max[1], node->box.max[2]);

            glEnd();

            // Set color for triangles
            glColor3f(1.0f, 0.0f, 0.0f);  // Red color

            glBegin(GL_TRIANGLES);
            for (const Triangle& triangle : node->triangles) {
                for (int i = 0; i < 3; ++i) {
                    glVertex3f(triangle.getVertices()[i][0], triangle.getVertices()[i][1], triangle.getVertices()[i][2]);
                }
            }
            glEnd();
        }

        // Recursively draw the bounding boxes of the child nodes
        if (node->left) {
            drawBoundingBoxesHelper(node->left);
        }
        if (node->right) {
            drawBoundingBoxesHelper(node->right);
        }
    }
    


};


#endif // KDTREE_H