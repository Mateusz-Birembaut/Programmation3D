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
    std::vector<unsigned int> indices;
};


class KdTree {

public:

    BoundingBox box;
    KdNode* root;
    int maxDepth;

    KdTree(const Mesh & mesh, int _maxDepth) {
        std::vector<unsigned int> indices(mesh.triangles.size());
        for (unsigned int i = 0; i < mesh.triangles.size(); ++i) {
            indices[i] = i;
        }
        maxDepth = _maxDepth;
        box = BoundingBox::meshBoundingBox(mesh, indices);
        root = new KdNode();
        root->node_box = box;
        KdTreeBuild( mesh ,box, root, indices);
    }

    void KdTreeBuild(const Mesh & _mesh ,BoundingBox& _parentBox, KdNode* _node, const std::vector<unsigned int>& _triangles_indexes, int depth = 0) {
        if (depth == maxDepth){ 
            _node->isLeaf = true;
            _node->indices = _triangles_indexes;
            _node->left = nullptr;
            _node->right = nullptr;
            return;
        }

        Vec3 diff = _parentBox.max - _parentBox.min;
        _node->dimSplit = diff.maxDimension(); 
        //_node->splitDistance = _parentBox.min[_node->dimSplit] + diff[_node->dimSplit] / 2; 
        _node->splitDistance = findBestSplit(_mesh ,_node, _parentBox, _triangles_indexes, _node->dimSplit, 10);

        std::vector<unsigned int> leftTriangles;
        std::vector<unsigned int> rightTriangles;

        for (const unsigned int & index : _triangles_indexes){ 
            const MeshTriangle& meshTriangle = _mesh.triangles[index];
            Triangle triangle(_mesh.vertices[meshTriangle[0]].position,
                            _mesh.vertices[meshTriangle[1]].position,
                            _mesh.vertices[meshTriangle[2]].position);
            BoundingBox triangleBox = BoundingBox::triangleBoundingBox(triangle);
            if (triangleBox.min[_node->dimSplit] <= _node->splitDistance) leftTriangles.push_back(index);
            if (triangleBox.max[_node->dimSplit] > _node->splitDistance) rightTriangles.push_back(index);

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
            KdTreeBuild( _mesh, leftBox, _node->left, leftTriangles, depth + 1);
            KdTreeBuild( _mesh, rightBox, _node->right, rightTriangles, depth + 1);
        }

    }

    float findBestSplit(const Mesh & _mesh,KdNode* _node, const BoundingBox& _parentBox, const std::vector<unsigned int >& _trianglesIndexes, int _dimSplit, int _splitsToTest) {
        Vec3 diff = _parentBox.max - _parentBox.min;
        float totalArea = BoundingBox::calculateSurfaceArea(_parentBox);

        float bestCost = FLT_MAX;
        float bestSplit = 0.0f;

        float step = diff[_dimSplit] / _splitsToTest;

        for (int i = 1; i < _splitsToTest; ++i) {
            float split = _parentBox.min[_dimSplit] + i * step;

            BoundingBox leftBox = _parentBox;
            leftBox.max[_dimSplit] = split;
            float leftArea = BoundingBox::calculateSurfaceArea(leftBox);

            BoundingBox rightBox = _parentBox;
            rightBox.min[_dimSplit] = split;
            float rightArea = BoundingBox::calculateSurfaceArea(rightBox);

            int numLeft = 0;
            int numRight = 0;


            for (const unsigned int & index : _trianglesIndexes) {
                const MeshTriangle& meshTriangle = _mesh.triangles[index];
                Triangle triangle(_mesh.vertices[meshTriangle[0]].position,
                                _mesh.vertices[meshTriangle[1]].position,
                                _mesh.vertices[meshTriangle[2]].position);
                BoundingBox triBox = BoundingBox::triangleBoundingBox(triangle);
                if (triBox.min[_dimSplit] <= split) numLeft ++;
                if (triBox.max[_dimSplit] > split) numRight ++;
            }

            float cost = calculateCost(numLeft, leftArea, numRight, rightArea, totalArea);

            if (cost < bestCost) {
                bestCost = cost;
                bestSplit = split;
            }
        }

        return bestSplit;
    }


    float calculateCost(int numLeft, float leftArea, int numRight, float rightArea, float totalArea) {
        return (numLeft * leftArea + numRight * rightArea) / totalArea;
    }


    RayTriangleIntersection traverse(const Mesh & _mesh, Ray const & ray, KdNode* _node, float t_start, float t_end) const {
        if (_node->isLeaf){
            RayTriangleIntersection intersection;
            for (const unsigned int & index : _node->indices) {
                const MeshTriangle& meshTriangle = _mesh.triangles[index];
                const Vec3& v0 = _mesh.vertices[meshTriangle[0]].position;
                const Vec3& v1 = _mesh.vertices[meshTriangle[1]].position;
                const Vec3& v2 = _mesh.vertices[meshTriangle[2]].position;
                RayTriangleIntersection tempIntersection = Triangle(v0, v1, v2).getIntersection(ray);
                if (tempIntersection.intersectionExists && tempIntersection.t < t_end && tempIntersection.t > t_start) {
                    intersection = tempIntersection;
                }
            }
            return intersection;
        }

        std::pair<float, float> interval = _node->node_box.intersect(ray);
        if (interval.first == INFINITY && interval.second == INFINITY) {
            return RayTriangleIntersection(); 
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
            return traverse(_mesh , ray, secondNode, t_start, t_end);
        } else if (t >= t_end) {
            return traverse(_mesh ,ray, firstNode, t_start, t_end);
        } else {
            RayTriangleIntersection hit = traverse(_mesh ,ray, firstNode, t_start, t);
            if (hit.intersectionExists && hit.t <= t && hit.t >= 0.0001) {
                return hit;
            }
            return traverse(_mesh ,ray, secondNode, t, t_end);
        }
    }


    void drawBoundingBoxesHelper(const KdNode * node) const {
        if (!node) return;
        if (node->left == nullptr && node->right == nullptr && node->indices.size() > 0) {

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