#ifndef OBJECTKDTREE_H
#define OBJECTKDTREE_H
#include <vecmath.h>

#include <map>
#include <vector>

#include "bound.hpp"
#include "hit.hpp"
#include "object3d.hpp"
using std::map;
using std::vector;

// Ref: ver.2020
class ObjectKDTreeNode {
public:
    Vector3f min, max;            // 坐标最小值和最大值
    vector<Object3D*>* faces;     // 子空间物体列表
    ObjectKDTreeNode *ls, *rs;    // 左右孩子
    // TODO: l, r 不知道什么用，待删
    int l, r;                     
    bool inside(Object3D* face) {
        Vector3f faceMin = face->min();
        Vector3f faceMax = face->max();
        // 范围位于包围盒内或者刚好贴在包围盒边缘，都算在包围盒内
        return (faceMin.x() < max.x() ||
                faceMin.x() == max.x() && faceMin.x() == faceMax.x()) &&
               (faceMax.x() > min.x() ||
                faceMax.x() == min.x() && faceMin.x() == faceMax.x()) &&
               (faceMin.y() < max.y() ||
                faceMin.y() == max.y() && faceMin.y() == faceMax.y()) &&
               (faceMax.y() > min.y() ||
                faceMax.y() == min.y() && faceMin.y() == faceMax.y()) &&
               (faceMin.z() < max.z() ||
                faceMin.z() == max.z() && faceMin.z() == faceMax.z()) &&
               (faceMax.z() > min.z() ||
                faceMax.z() == min.z() && faceMin.z() == faceMax.z());
    }
};

// Ref: ver.2020
// TODO：修改建树
class ObjectKDTree {
    int n;
    Vector3f** vertices;  // 顶点数组
    // 按照 KD-tree 算法建树
    ObjectKDTreeNode* build(int depth, int axis, vector<Object3D*>* faces,
                            const Vector3f& min, const Vector3f& max) {
        ObjectKDTreeNode* p = new ObjectKDTreeNode;
        p->min = min;
        p->max = max;
        Vector3f maxL, minR;  // 左子树最小坐标与右子树最大坐标
        if (axis == 0) {
            maxL =
                Vector3f((p->min.x() + p->max.x()) / 2, p->max.y(), p->max.z());
            minR =
                Vector3f((p->min.x() + p->max.x()) / 2, p->min.y(), p->min.z());
        } else if (axis == 1) {
            maxL =
                Vector3f(p->max.x(), (p->min.y() + p->max.y()) / 2, p->max.z());
            minR =
                Vector3f(p->min.x(), (p->min.y() + p->max.y()) / 2, p->min.z());
        } else {
            maxL =
                Vector3f(p->max.x(), p->max.y(), (p->min.z() + p->max.z()) / 2);
            minR =
                Vector3f(p->min.x(), p->min.y(), (p->min.z() + p->max.z()) / 2);
        }
        p->faces = new vector<Object3D*>;
        for (auto face : *faces)
            if (p->inside(face)) p->faces->push_back(face);

        const int max_faces = 128;  // 叶节点最大容量
        const int max_depth = 24;   // 最大树深

        // 递归建树
        if (p->faces->size() > max_faces && depth < max_depth) {
            p->ls = build(depth + 1, (axis + 1) % 3, p->faces, min, maxL);
            p->rs = build(depth + 1, (axis + 1) % 3, p->faces, minR, max);

            // 横跨分界面的节点归给父亲    
            vector<Object3D*>*faceL = p->ls->faces, *faceR = p->rs->faces;
            map<Object3D*, int> cnt;
            for (auto face : *faceL) cnt[face]++;
            for (auto face : *faceR) cnt[face]++;
            p->ls->faces = new vector<Object3D*>;
            p->rs->faces = new vector<Object3D*>;
            p->faces->clear();
            for (auto face : *faceL)
                if (cnt[face] == 1)
                    p->ls->faces->push_back(face);
                else
                    p->faces->push_back(face);
            for (auto face : *faceR)
                if (cnt[face] == 1) p->rs->faces->push_back(face);
        } else
            p->ls = p->rs = nullptr;
        return p;
    }

    // 递归获取所有物体
    void getFaces(ObjectKDTreeNode* p, vector<Object3D*>* faces) {
        for (auto face : *(p->faces)) faces->push_back(face);
        // TODO: l, r 不知道什么用，待删
        p->l = faces->size();
        p->r = faces->size();
        if (p->ls) getFaces(p->ls, faces);
        if (p->rs) getFaces(p->rs, faces);
    }

public:
    ObjectKDTreeNode* root;    // 根节点
    vector<Object3D*>* faces;  // 所有物体列表
    ObjectKDTree(vector<Object3D*>* faces) {
        Vector3f min = Vector3f(INF, INF, INF);
        Vector3f max = -min;
        // 对每个分量取最小最大值作为边界
        for (auto face : *faces) {
            min = minE(min, face->min());
            max = maxE(max, face->max());
        }
        // 建树（初始深度为 1，轴为 x 轴）
        root = build(1, 0, faces, min, max);
        std::cout << "建树完成" << std::endl;
        // TODO：为什么这里 this->faces 不直接等于 faces 呢？
        this->faces = new vector<Object3D*>;
        getFaces(root, this->faces);    
    }

    // 计算光线与给定子空间的交点距离（注意：这是和包围盒相交，不见得和包围盒内物体相交）
    float cuboidIntersect(ObjectKDTreeNode* p, const Ray& ray) const {
        float t = INF;
        if (!p) return t;
        AABB(p->min, p->max).intersect(ray, t);
        return t;
    }

    // 检测光线与整个场景是否有交点，并更新最近的交点信息至 hit 中
    bool intersect(const Ray& ray, Hit& hit, float tmin = TMIN) const {
        Object3D* nextFace = nullptr;
        return intersect2(root, ray, nextFace, hit);
    }

    // 递归地在 KD 树中查找射线与物体的交点
    bool intersect2(ObjectKDTreeNode* p, const Ray& ray, Object3D*& nextFace,
                   Hit& hit) const {
        bool flag = false;  // 表示有无交点
        for (int i = 0; i < p->faces->size(); ++i)
            if ((*p->faces)[i]->intersect(ray, hit, TMIN)) {
                nextFace = (*p->faces)[i];
                flag = true;
            }
        float tl = cuboidIntersect(p->ls, ray),
              tr = cuboidIntersect(p->rs, ray);
        // 根据和包围盒交点深浅决定先递归检查哪个
        if (tl < tr) {
            if (hit.getT() <= tl) return flag;
            if (p->ls) flag |= intersect2(p->ls, ray, nextFace, hit);
            if (hit.getT() <= tr) return flag;
            if (p->rs) flag |= intersect2(p->rs, ray, nextFace, hit);
        } else {
            if (hit.getT() <= tr) return flag;
            if (p->rs) flag |= intersect2(p->rs, ray, nextFace, hit);
            if (hit.getT() <= tl) return flag;
            if (p->ls) flag |= intersect2(p->ls, ray, nextFace, hit);
        }
        return flag;
    }
};

#endif  // !OBJECTKDTREE_H