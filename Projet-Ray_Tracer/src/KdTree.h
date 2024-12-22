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

    KdTreeNode() : left(nullptr), right(nullptr) {}

    KdTreeNode* KdTreeBuild(const std::vector<Triangle>& _triangles, int _depth){
        if (_triangles.size() == 0){
            std::cerr << "Error: No triangles for this node" << std::endl;
            return nullptr;
        }

        if (_depth == 0 || _triangles.size() <= 2) {
            KdTreeNode* node = new KdTreeNode();
            node->isLeaf = true;
            for (int i = 0; i < _triangles.size(); i++){
                node->triangles.push_back(_triangles[i]);
            }
            std::cout << "Leaf : " << node->triangles.size() << std::endl;
            return node;

        }else {


            KdTreeNode* node = new KdTreeNode();

            BoundingBox box = BoundingBox().meshBoundingBox(_triangles);

            Vec3 diff = box.max - box.min;
            node->dimSplit = diff.maxDimension();

            node->splitDistance = box.min[node->dimSplit] + diff[node->dimSplit] / 2; // pour l'instant coupe au milieu

            std::vector<Triangle> leftTriangles;
            std::vector<Triangle> rightTriangles;

            for (int i = 0; i < _triangles.size(); i++){ 
                Triangle triangle = _triangles[i];
                BoundingBox triangleBox = BoundingBox().triangleBoundingBox(triangle);

                
                if (triangleBox.max[node->dimSplit] < node->splitDistance){
                    leftTriangles.push_back(triangle);
                } else if (triangleBox.min[node->dimSplit] > node->splitDistance){
                    rightTriangles.push_back(triangle);
                } else { // quand le triangle est coupé en deux on le met dans les deux
                    leftTriangles.push_back(triangle);
                    rightTriangles.push_back(triangle);
                }
                
            }

            std::cout << "Left: " << leftTriangles.size() << " Right: " << rightTriangles.size() << std::endl;

            node->left = KdTreeBuild(leftTriangles, _depth - 1);
            node->right = KdTreeBuild(rightTriangles, _depth - 1);

            return node;
        }
    }

    RayTriangleIntersection traverse(Ray const & ray, RayTriangleIntersection & intersection) const { 
        //TODO : Implementer la fonction de traversée de l'arbre
    }


    void testKdTreeBuild() {
        std::vector<Triangle> triangles;
        triangles.push_back(Triangle(Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0)));
        triangles.push_back(Triangle(Vec3(1, 1, 1), Vec3(2, 1, 1), Vec3(1, 2, 1)));
        triangles.push_back(Triangle(Vec3(2, 2, 2), Vec3(3, 2, 2), Vec3(2, 3, 2)));
        //triangles.push_back(Triangle(Vec3(12, 12, 12), Vec3(11, 11, 11), Vec3(10, 10, 10)));


        KdTreeNode* root = new KdTreeNode();
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

        KdTreeNode* root = new KdTreeNode();
        root = root->KdTreeBuild(testTriangles, 3);

        Ray ray(Vec3(0, 0, -1), Vec3(0, 0, 1));
        RayTriangleIntersection intersection;
        intersection.intersectionExists = false;
        intersection.t = std::numeric_limits<float>::max();

        root->traverse(ray, intersection);
        if (intersection.intersectionExists)
        {
            std::cout << "Intersection found: " << intersection.intersectionExists 
                << " at distance " << intersection.t << std::endl;
        }
    
    }   


    void testKdTreeBuildComplex() {
        std::vector<Triangle> triangles;

        triangles.push_back(Triangle(Vec3(-1, -1, 0), Vec3(-0.5, -1, 0), Vec3(-1, -0.5, 0)));
        triangles.push_back(Triangle(Vec3(0, 0, 0),   Vec3(1, 0, 0),    Vec3(0, 1, 0)));
        triangles.push_back(Triangle(Vec3(2, 2, 0),   Vec3(3, 2, 0),    Vec3(2, 3, 0)));
        triangles.push_back(Triangle(Vec3(5, 5, 5),   Vec3(4, 6, 5),    Vec3(5, 6, 5)));
        triangles.push_back(Triangle(Vec3(10, 10, -2),Vec3(12, 10, -2), Vec3(10, 12, -2)));
        triangles.push_back(Triangle(Vec3(-2, 2, 2),  Vec3(-3, 2, 2),   Vec3(-2, 3, 2)));

        KdTreeNode* root = new KdTreeNode();

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
        }
    } 


};


#endif // KDTREE_H