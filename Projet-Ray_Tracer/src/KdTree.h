#ifndef KDTREE_H
#define KDTREE_H

#include "BoundingBox.h"
#include <vector>
#include "BoundingBox.h"
#include "Sphere.h"
#include "Square.h"

struct KdNode {
    KdNode *left, *right;

    BoundingBox node_box;
    int dimSplit;
    float splitDistance;

    bool isLeaf = false;
    std::vector<Triangle> triangles;
};


class KdTree {


public:

    BoundingBox box;
    KdNode* root;
    int maxDepth;

    KdTree(const Mesh * mesh, int _maxDepth) {
        std::vector<Triangle> triangles;
        for (const auto& meshTriangle : mesh->triangles) {
            Vec3 v0 = mesh->vertices[meshTriangle[0]].position;
            Vec3 v1 = mesh->vertices[meshTriangle[1]].position;
            Vec3 v2 = mesh->vertices[meshTriangle[2]].position;
            triangles.emplace_back(v0, v1, v2);
        }
        std::cout << "Building KdTree with " << triangles.size() << " triangles" << std::endl;
        maxDepth = _maxDepth;
        box = BoundingBox().meshBoundingBox(triangles);
        root = new KdNode();
        root->node_box = box;
        KdTreeBuild(box, root, triangles);
    }

    void KdTreeBuild(BoundingBox _parentBox, KdNode* _node, const std::vector<Triangle>& _triangles, int depth = 0) {
        if (depth == maxDepth || _triangles.size() == 0){ // depth max atteint
            _node->isLeaf = true;
            _node->triangles = _triangles;
            _node->left = nullptr;
            _node->right = nullptr;
            return;
        }

        //std::cout << "Bounding box: " << _parentBox.min << " " << _parentBox.max << std::endl;
        
        Vec3 diff = _parentBox.max - _parentBox.min;
        _node->dimSplit = diff.maxDimension(); // quelle dimension a split
        _node->splitDistance = _parentBox.min[_node->dimSplit] + diff[_node->dimSplit] / 2; // on coupe en deux

        std::vector<Triangle> leftTriangles;
        std::vector<Triangle> rightTriangles;

        for (const Triangle& triangle : _triangles){ 
            BoundingBox triangleBox = BoundingBox().triangleBoundingBox(triangle);
            if (triangleBox.max[_node->dimSplit] < _node->splitDistance){
                leftTriangles.push_back(triangle);
            } else if (triangleBox.min[_node->dimSplit] > _node->splitDistance){
                rightTriangles.push_back(triangle);
            } else {
                leftTriangles.push_back(triangle);
                rightTriangles.push_back(triangle);
            }
        }
        _node->left = new KdNode();
        _node->right = new KdNode();

        BoundingBox leftBox = _parentBox;
        leftBox.max[_node->dimSplit] = _node->splitDistance;
        BoundingBox rightBox = _parentBox;
        rightBox.min[_node->dimSplit] = _node->splitDistance;

        _node->left->node_box = leftBox;
        _node->right->node_box = rightBox;

        if (depth != maxDepth){
            KdTreeBuild(leftBox, _node->left, leftTriangles, depth + 1);
            KdTreeBuild(rightBox, _node->right, rightTriangles, depth + 1);
        }

    }

    RayTriangleIntersection traverse(Ray const & ray, KdNode* _node, float t_start, float t_end) const {
        if (_node->isLeaf){
            RayTriangleIntersection intersection;
            intersection.intersectionExists = false;
            intersection.t = FLT_MAX;
            for (const Triangle& triangle : _node->triangles) {
                RayTriangleIntersection tempIntersection = triangle.getIntersection(ray);
                if (tempIntersection.intersectionExists && tempIntersection.t < t_end && tempIntersection.t > t_start) {
                    intersection = tempIntersection;
                }
            }
            return intersection;
        }

        float t = (_node->splitDistance - ray.origin()[_node->dimSplit]) / ray.direction()[_node->dimSplit];

        KdNode* firstNode;
        KdNode* secondNode;

        if (ray.direction()[_node->dimSplit] >= 0) {
            firstNode = _node->left;
            secondNode = _node->right;
        } else {
            firstNode = _node->right;
            secondNode = _node->left;
        }

        if (t <= t_start) {
            return traverse(ray, secondNode, t_start, t_end);
        } else if (t >= t_end) {
            return traverse(ray, firstNode, t_start, t_end);
        } else {
            RayTriangleIntersection hit = traverse(ray, firstNode, t_start, t);
            if (hit.intersectionExists && hit.t <= t && hit.t >= 0.0001) {
                return hit;
            }
            return traverse(ray, secondNode, t, t_end);
        }
    }




    void drawBoundingBoxesHelper(const KdNode* node) const {
        if (!node) return;
        if (node->left == nullptr && node->right == nullptr && node->triangles.size() > 0) {

            glColor3f(1.0f, 1.0f, 1.0f);   

            glLineWidth(1.0f);

            glBegin(GL_LINES);

            glVertex3f(node->node_box.min[0], node->node_box.min[1], node->node_box.min[2]);
            glVertex3f(node->node_box.max[0], node->node_box.min[1], node->node_box.min[2]);

            glVertex3f(node->node_box.min[0], node->node_box.min[1], node->node_box.min[2]);
            glVertex3f(node->node_box.min[0], node->node_box.max[1], node->node_box.min[2]);

            glVertex3f(node->node_box.min[0], node->node_box.min[1], node->node_box.min[2]);
            glVertex3f(node->node_box.min[0], node->node_box.min[1], node->node_box.max[2]);

            glVertex3f(node->node_box.max[0], node->node_box.max[1], node->node_box.max[2]);
            glVertex3f(node->node_box.min[0], node->node_box.max[1], node->node_box.max[2]);

            glVertex3f(node->node_box.max[0], node->node_box.max[1], node->node_box.max[2]);
            glVertex3f(node->node_box.max[0], node->node_box.min[1], node->node_box.max[2]);

            glVertex3f(node->node_box.max[0], node->node_box.max[1], node->node_box.max[2]);
            glVertex3f(node->node_box.max[0], node->node_box.max[1], node->node_box.min[2]);

            glVertex3f(node->node_box.min[0], node->node_box.max[1], node->node_box.min[2]);
            glVertex3f(node->node_box.max[0], node->node_box.max[1], node->node_box.min[2]);

            glVertex3f(node->node_box.min[0], node->node_box.max[1], node->node_box.min[2]);
            glVertex3f(node->node_box.min[0], node->node_box.max[1], node->node_box.max[2]);

            glVertex3f(node->node_box.max[0], node->node_box.min[1], node->node_box.min[2]);
            glVertex3f(node->node_box.max[0], node->node_box.max[1], node->node_box.min[2]);

            glVertex3f(node->node_box.max[0], node->node_box.min[1], node->node_box.min[2]);
            glVertex3f(node->node_box.max[0], node->node_box.min[1], node->node_box.max[2]);

            glVertex3f(node->node_box.min[0], node->node_box.min[1], node->node_box.max[2]);
            glVertex3f(node->node_box.max[0], node->node_box.min[1], node->node_box.max[2]);

            glVertex3f(node->node_box.min[0], node->node_box.min[1], node->node_box.max[2]);
            glVertex3f(node->node_box.min[0], node->node_box.max[1], node->node_box.max[2]);

            glEnd();

            glColor3f(1.0f, 0.0f, 0.0f); 

            glBegin(GL_TRIANGLES);
            for (const Triangle& triangle : node->triangles) {
                for (int i = 0; i < 3; ++i) {
                    glVertex3f(triangle.getVertices()[i][0], triangle.getVertices()[i][1], triangle.getVertices()[i][2]);
                }
            }
            glEnd();
        }

        if (node->left) {
            drawBoundingBoxesHelper(node->left);
        }
        if (node->right) {
            drawBoundingBoxesHelper(node->right);
        }
    }

};



#endif