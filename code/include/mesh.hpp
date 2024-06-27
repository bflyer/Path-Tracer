#ifndef MESH_H
#define MESH_H

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

#include "object3d.hpp"
#include "triangle.hpp"
#include "Vector2f.h"
#include "Vector3f.h"
#include "obj_kdtree.hpp"
#include "bound.hpp"
#include "ray.hpp"
#include "utils.hpp"

// Ref: ver.2020
class Mesh : public Object3D {

public:
    // TODO: 简化构造函数
    Mesh(const char *filename, Material *m) : Object3D(m) {
        std::ifstream f;                    
        f.open(filename);
        if (!f.is_open()) {
            std::cout << "Cannot open " << filename << "\n";
            return;
        }
        std::vector<TriangleIndex> vIdx, tIdx, nIdx;  // 面的顶点索引
        std::vector<Vector3f> v, vn;         // 顶点集，法线索引集
        std::vector<Vector2f> vt;            // 纹理坐标集
        
        std::string line;
        std::string bslash("/"), space(" ");
        std::string token;
        int texID;
        while (true) {
            std::getline(f, line);
            // 跳过空行、注释行（以#开头）和长度小于3的行。
            if (line.size() < 3 || line.at(0) == '#')
                continue;

            std::stringstream ss(line);
            ss >> token;
            // 若关键词为 v，则读取顶点坐标并更新 AABB
            if (token == "v") {
                Vector3f vec;
                ss >> vec[0] >> vec[1] >> vec[2];
                v.push_back(vec);
                aabb.updateBound(vec);
            }
            // 若关键词为 f，则读取面顶点索引
            else if (token == "f") {
                // 存储顶点索引、纹理坐标索引和法线索引
                TriangleIndex vId, tId, nId;
                for (int i = 0; i < 3; ++i) {
                    std::string str;
                    ss >> str;
                    std::vector<std::string> id = split(str, bslash);
                    vId[i] = atoi(id[0].c_str()) - 1;
                    if (id.size() > 1) {
                        tId[i] = atoi(id[1].c_str()) - 1;
                    }
                    if (id.size() > 2) {
                        nId[i] = atoi(id[2].c_str()) - 1;
                    }
                }
                vIdx.push_back(vId);
                tIdx.push_back(tId);
                nIdx.push_back(nId);
            } 
            // 若关键词为 vt，则读取纹理坐标
            else if (token == "vt") {
                Vector2f texcoord;
                ss >> texcoord[0];
                ss >> texcoord[1];
                vt.push_back(texcoord);
            } 
            // 若关键词为 vn，则读取法线坐标
            else if (token == "vn") {
                Vector3f vec;
                ss >> vec[0] >> vec[1] >> vec[2];
                vn.push_back(vec);
            }
            if (f.eof()) {
                break;
            }
        }
        f.close();
        for (int triId = 0; triId < (int)vIdx.size(); ++triId) {
            // A. 创建一个三角形并入队
            TriangleIndex &vIndex = vIdx[triId];
            triangles.push_back((Object3D *)new Triangle(
                v[vIndex[0]], v[vIndex[1]], v[vIndex[2]], m));

            // B. 若纹理坐标合规，则设置纹理坐标
            TriangleIndex &tIndex = tIdx[triId];
            if (tIndex.valid())
                ((Triangle *)triangles.back())
                    ->setVT(vt[tIndex[0]], vt[tIndex[1]], vt[tIndex[2]]);

            // C. 若法线合规，则设置法线
            TriangleIndex &nIndex = nIdx[triId];
            if (nIndex.valid())
                ((Triangle *)triangles.back())
                    ->setVNorm(vn[nIndex[0]], vn[nIndex[1]], vn[nIndex[2]]);
            // }
        }
        area = 0;
        for (auto triangle : triangles) {
            area += triangle->getArea();
        }
        kdTree = new ObjectKDTree(&triangles);
    }

    ~Mesh() {
        for (auto triangle : triangles) delete triangle;
        delete kdTree;
    }

    vector<string> split(const string& str, const string& pattern) {
        vector<string> result;
        std::string::size_type prev = 0, pos;

        // 搜寻分隔符直至找不到
        while ((pos = str.find(pattern, prev)) != std::string::npos) {
            // 将 prev 到 pos 之前的子串添加到结果中
            result.push_back(str.substr(prev, pos - prev));
            // 更新 prev 到分隔符之后的位置
            prev = pos + pattern.size();
        }

        // 添加最后一个子串
        if (prev < str.size()) {
            result.push_back(str.substr(prev));
        }

        return result;
    }

    struct TriangleIndex {
        TriangleIndex() {
            x[0] = -1; x[1] = -1; x[2] = -1;
        }
        int &operator[](const int i) { return x[i]; }
        // By Computer Graphics convention, counterclockwise winding is front face
        int x[3]{};
        bool valid() { return x[0] != -1 && x[1] != -1 && x[2] != -1; }
    };

    std::vector<Object3D *> triangles;
    // KD-tree 与 包围盒加速求交
    bool intersect(const Ray &ray, Hit &hit, double tmin) override {
        double t;  // 与包围盒相交的 t
        // 若与包围盒无交，直接返回 false
        if (!aabb.intersect(ray, t)) return false;
        // 若与包围盒的交点不如当前交点近，也返回 false
        if (t > hit.getT()) return false;
        // KD-tree 加速求交
        return kdTree->intersect(ray, hit);
    }

    // // 顺序遍历求交
    // bool intersect (const Ray &ray, Hit &hit, double tmin) override {
    //     bool flag = false;
    //     for (auto triangle : triangles) flag |= triangle->intersect(ray, hit, TMIN);
    //     return flag;
    // }

    // 顺序遍历求交
    bool sequentialSearch(const Ray &ray, Hit &hit) {
        bool flag = false;
        for (auto triangle : triangles) flag |= triangle->intersect(ray, hit, TMIN);
        return flag;
    }

    Vector3f sample() const override {
        int id = (int)(RAND2 * triangles.size());
        return triangles[id]->sample();
    }

    Vector3f min() const override { return aabb.bounds[0]; }
    Vector3f max() const override { return aabb.bounds[1]; }
    Vector3f center() const override {
        return (aabb.bounds[0] + aabb.bounds[1]) / 2;
    }
    vector<Object3D *> getFaces() override { return {(Object3D *)this}; }
    double getArea() const override { return area; }

    // TODO: 待删
    std::vector<Vector3f> v;
    std::vector<TriangleIndex> t;
    std::vector<Vector3f> n;
    

private:
    AABB aabb;
    ObjectKDTree *kdTree;
    double area;
    // Normal can be used for light estimation
    void computeNormal();
};

#endif
