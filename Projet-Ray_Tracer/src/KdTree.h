#ifndef KDTREE_H
#define KDTREE_H

#include "BoundingBox.h"
#include <vector>
#include "BoundingBox.h"
#include "Sphere.h"
#include "Square.h"
#include "Photon.h"
#include <queue>
#include <chrono>

template <typename T>
struct KdNode {
    KdNode *left, *right;

    BoundingBox node_box;
    int dimSplit;
    float splitDistance;

    bool isLeaf = false;
    std::vector<T> primitives;
};


class KdTree {


public:

    BoundingBox box;
    KdNode<Triangle>* root;
    int maxDepth;

    KdTree(const Mesh * mesh, int _maxDepth) {
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<Triangle> triangles;
        for (const auto& meshTriangle : mesh->triangles) {
            Vec3 v0 = mesh->vertices[meshTriangle[0]].position;
            Vec3 v1 = mesh->vertices[meshTriangle[1]].position;
            Vec3 v2 = mesh->vertices[meshTriangle[2]].position;
            triangles.emplace_back(v0, v1, v2);
        }
        std::cout << "Building KdTree with " << triangles.size() << " triangles" << std::endl;
        maxDepth = _maxDepth;
        box = BoundingBox::meshBoundingBox(triangles);
        root = new KdNode<Triangle>();
        root->node_box = box;
        KdTreeBuild(box, root, triangles);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "KdTree built in " << elapsed.count() << " seconds. \n" << std::endl;
    }

    template <typename T>
    void KdTreeBuild(BoundingBox _parentBox, KdNode<T>* _node, const std::vector<Triangle>& _triangles, int depth = 0) {
        if (depth == maxDepth || _triangles.size() <= 1) { 
            _node->isLeaf = true;
            _node->primitives = _triangles;
            _node->left = nullptr;
            _node->right = nullptr;
            return;
        }

        Vec3 diff = _parentBox.max - _parentBox.min;
        _node->dimSplit = diff.maxDimension(); 
        //_node->splitDistance = _parentBox.min[_node->dimSplit] + diff[_node->dimSplit] / 2; 
        _node->splitDistance = findBestSplit(_node, _parentBox, _triangles, _node->dimSplit, 10);

        std::vector<Triangle> leftTriangles;
        std::vector<Triangle> rightTriangles;

        for (const Triangle& triangle : _triangles){ 
            BoundingBox triangleBox = BoundingBox::triangleBoundingBox(triangle);
            if (triangleBox.min[_node->dimSplit] <= _node->splitDistance) leftTriangles.push_back(triangle);
            if (triangleBox.max[_node->dimSplit] > _node->splitDistance) rightTriangles.push_back(triangle);

        }
        _node->left = new KdNode<Triangle>();
        _node->right = new KdNode<Triangle>();

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

    template <typename T>
    float findBestSplit(KdNode<T>* _node, const BoundingBox& _parentBox, const std::vector<Triangle>& _triangles, int _dimSplit, int _splitsToTest) {
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


            for (const Triangle& triangle : _triangles) {
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


    template <typename T>
    RayTriangleIntersection traverse(Ray const & ray, KdNode<T>* _node, float t_start, float t_end) const {
        if (_node->isLeaf){
            RayTriangleIntersection intersection;
            for (const Triangle& triangle : _node->primitives) {
                RayTriangleIntersection tempIntersection = triangle.getIntersection(ray);
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

        KdNode<Triangle>* firstNode;
        KdNode<Triangle>* secondNode;

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


    template <typename T>
    void drawBoundingBoxesHelper(const KdNode<T>* node, bool showBoxes) const {
        if (!node) return;
        if (node->left == nullptr && node->right == nullptr && node->triangles.size() > 0) {

            if (showBoxes) {
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
            }

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
            drawBoundingBoxesHelper(node->left, showBoxes);
        }
        if (node->right) {
            drawBoundingBoxesHelper(node->right, showBoxes);
        }
    }

};

class KdTreePhotonMap {

public:

    BoundingBox box;
    KdNode<Photon>* root;
    int maxDepth;

    KdTreePhotonMap(std::vector<Photon> & photons, int _maxDepth) {
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "KdTree Photon map with " << photons.size() << " photons" << std::endl;
        maxDepth = _maxDepth;
        box = BoundingBox::photonMapBoundingBox(photons);
        root = new KdNode<Photon>();
        root->node_box = box;
        KdTreePhotonMapBuild(box, root, photons);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "KdTree Photon map built in " << elapsed.count() << " seconds. \n" << std::endl;
    }

    void KdTreePhotonMapBuild(BoundingBox _parentBox, KdNode<Photon>* _node, std::vector<Photon>& _photons, int depth = 0){
        if (depth == maxDepth || _photons.size() <= 1){ 
            _node->isLeaf = true;
            _node->primitives = _photons;
            _node->left = nullptr;
            _node->right = nullptr;
            return;
        }

        Vec3 diff = _parentBox.max - _parentBox.min;
        _node->dimSplit = diff.maxDimension(); 
        //_node->splitDistance = _parentBox.min[_node->dimSplit] + diff[_node->dimSplit] / 2; 
        _node->splitDistance = findMedian(_photons, _node->dimSplit);

        std::vector<Photon> leftPhotons;
        std::vector<Photon> rightPhotons;

        for (const Photon & photon : _photons){ 
            if (photon.position[_node->dimSplit] <= _node->splitDistance) leftPhotons.push_back(photon);
            if (photon.position[_node->dimSplit] > _node->splitDistance) rightPhotons.push_back(photon);

        }
        
        _node->left = new KdNode<Photon>();
        _node->right = new KdNode<Photon>();

        BoundingBox leftBox = _parentBox;
        leftBox.max[_node->dimSplit] = _node->splitDistance;
        BoundingBox rightBox = _parentBox;
        rightBox.min[_node->dimSplit] = _node->splitDistance;

        _node->left->node_box = leftBox;
        _node->right->node_box = rightBox;

        if (depth != maxDepth){
            KdTreePhotonMapBuild(leftBox, _node->left, leftPhotons, depth + 1);
            KdTreePhotonMapBuild(rightBox, _node->right, rightPhotons, depth + 1);
        }

    }

    std::vector<Photon> findNearestPhotons(const Vec3& position, float radius) {
        std::vector<Photon> result;

        BoundingBox cercleBox = BoundingBox::createSphereBox(position, radius);

        findNearestPhotonsRecursive(root, cercleBox , position, radius, result);

        return result;
    }

    void findNearestPhotonsRecursive(KdNode<Photon>* _node,const BoundingBox & cercleBox, const Vec3& _position, float _radius, std::vector<Photon>& _nearestPhotons) {
        if (!_node) return;
        if (_node->isLeaf){
            for (const Photon & photon : _node->primitives) {
                float distance = (photon.position - _position).length();
                if (distance < _radius) {
                    _nearestPhotons.push_back(photon);
                }
            }
            return;
        }

        if (_node->left->node_box.overlaps(cercleBox)){ 
            findNearestPhotonsRecursive(_node->left, cercleBox, _position, _radius, _nearestPhotons);
        }

        if (_node->right->node_box.overlaps(cercleBox)){
            findNearestPhotonsRecursive(_node->right, cercleBox, _position, _radius, _nearestPhotons);
        }
    }


    float findMedian(std::vector<Photon>& photons, int dimSplit) {
        size_t n = photons.size();
        if (n == 0) return 0.0f;

        std::nth_element(photons.begin(), photons.begin() + n / 2, photons.end(),
                        [dimSplit](const Photon& a, const Photon& b) {
                            return a.position[dimSplit] < b.position[dimSplit];
                        });

        return photons[n / 2].position[dimSplit];
    }

};




/*     void findNearestPhotonsRecursive(KdNode<Photon>* _node, const Vec3& _position, unsigned int _k, std::priority_queue<PhotonDistance>& _nearestPhotons) {
        if (!_node) return;

        if (_node->isLeaf) {
            for (const Photon & photon : _node->primitives) {
                float distance = (photon.position - _position).length();
                if (_nearestPhotons.size() < _k) {
                    _nearestPhotons.push(PhotonDistance(photon, distance));
                } else {
                    if (distance < _nearestPhotons.top().distance) {
                        _nearestPhotons.pop();
                        _nearestPhotons.push(PhotonDistance(photon, distance));
                    }
                }
            }
            return;
        }

        int splitDim = _node->dimSplit;
        KdNode<Photon>* nearNode = _node->left;
        KdNode<Photon>* farNode = _node->right;

        if (_position[splitDim] > _node->splitDistance) {
            findNearestPhotonsRecursive(nearNode, _position, _k, _nearestPhotons);
        }

        if (farNode && (_nearestPhotons.size() < _k || std::abs(_position[splitDim] - _node->splitDistance) < _nearestPhotons.top().distance)) {
            findNearestPhotonsRecursive(farNode, _position, _k, _nearestPhotons);
        }
    } */




#endif