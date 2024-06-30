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

#define less(fMin1, bMax, fMax1) (fMin1 < bMax || (fMin1 == bMax && fMin1 == fMax1))
#define greater(fMax2, bMin, fMin2) (fMax2 > bMin || (fMax2 == bMin && fMax2 == fMin2))

// Ref: ver.2020
class ObjectKDTreeNode {
public:
    Vector3f min, max;            // 坐标最小值和最大值
    vector<Object3D*>* faces;     // 子空间物体列表
    ObjectKDTreeNode *l, *r;    // 左右孩子

    bool inside(Object3D* face) {
        Vector3f faceMin = face->min();
        Vector3f faceMax = face->max();
        // 范围位于包围盒内、刚好贴在包围盒边缘或横跨包围盒，都算在包围盒内
        return less(faceMin.x(), max.x(), faceMax.x()) && greater(faceMax.x(), min.x(), faceMin.x()) && 
               less(faceMin.y(), max.y(), faceMax.y()) && greater(faceMax.y(), min.y(), faceMin.y()) && 
               less(faceMin.z(), max.z(), faceMax.z()) && greater(faceMax.z(), min.z(), faceMin.z());
    }
};

// Ref: ver.2020
class ObjectKDTree {
    int n;
    Vector3f** vertices;  // 顶点数组
    // 按照 KD-tree 算法建树
    ObjectKDTreeNode* build(int depth, int axis, vector<Object3D*>* faces,
                            const Vector3f& min, const Vector3f& max) {
        ObjectKDTreeNode* p = new ObjectKDTreeNode;
        p->min = min;
        p->max = max;
        Vector3f maxL, minR;  // 左子树最大坐标与右子树最小坐标
        // 按轴选分界面
        if (axis == 0) {
            maxL = Vector3f((p->min.x() + p->max.x()) / 2, p->max.y(), p->max.z());
            minR = Vector3f((p->min.x() + p->max.x()) / 2, p->min.y(), p->min.z());
        } else if (axis == 1) {
            maxL = Vector3f(p->max.x(), (p->min.y() + p->max.y()) / 2, p->max.z());
            minR = Vector3f(p->min.x(), (p->min.y() + p->max.y()) / 2, p->min.z());
        } else {
            maxL = Vector3f(p->max.x(), p->max.y(), (p->min.z() + p->max.z()) / 2);
            minR = Vector3f(p->min.x(), p->min.y(), (p->min.z() + p->max.z()) / 2);
        }
        p->faces = new vector<Object3D*>;
        // 将所有属于该空间的面片放入该空间的面片集合
        for (auto face : *faces)
            if (p->inside(face)) p->faces->push_back(face);

        const int max_faces = 64;  // 叶节点最大容量
        const int max_depth = 26;   // 最大树深

        // 递归建树
        if (p->faces->size() > max_faces && depth < max_depth) {
            // 按顺序选择下一个轴交给子树建树
            p->l = build(depth + 1, (axis + 1) % 3, p->faces, min, maxL);
            p->r = build(depth + 1, (axis + 1) % 3, p->faces, minR, max);

            // 横跨分界面的节点归给父亲    
            vector<Object3D*>*faceL = p->l->faces, *faceR = p->r->faces;
            // 记录每个面片在左右子树中出现的总次数
            map<Object3D*, int> cnt;
            for (auto face : *faceL) cnt[face]++;
            for (auto face : *faceR) cnt[face]++;
            p->l->faces = new vector<Object3D*>;
            p->r->faces = new vector<Object3D*>;
            p->faces->clear();

            // 重新分配面片，使其在任何一个空间中只可能出现一次
            for (auto face : *faceL)
                // 只出现一次的留在子空间
                if (cnt[face] == 1)
                    p->l->faces->push_back(face);
                // 出现了两次的（即横跨分界面）交给父亲
                else
                    p->faces->push_back(face);
            for (auto face : *faceR)
                if (cnt[face] == 1) p->r->faces->push_back(face);
        } else
            p->l = p->r = nullptr;
        return p;
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
        this->faces = faces;
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
        // 1. 先在本空间内进行求交检测
        for (int i = 0; i < p->faces->size(); ++i)
            if ((*p->faces)[i]->intersect(ray, hit, TMIN)) {
                nextFace = (*p->faces)[i];
                flag = true;
            }
        // 2. 再递归地在左右子空间中进行求交检测
        float tl = cuboidIntersect(p->l, ray),
              tr = cuboidIntersect(p->r, ray);
        // 根据射线和包围盒交点的深浅决定先递归检查哪个
        if (tl < tr) {
            // 3. 若当前最近交点比射线与包围盒的交点还近，直接返回
            if (hit.getT() <= tl) return flag;
            if (p->l) flag |= intersect2(p->l, ray, nextFace, hit);
            if (hit.getT() <= tr) return flag;
            if (p->r) flag |= intersect2(p->r, ray, nextFace, hit);
        } else {
            if (hit.getT() <= tr) return flag;
            if (p->r) flag |= intersect2(p->r, ray, nextFace, hit);
            if (hit.getT() <= tl) return flag;
            if (p->l) flag |= intersect2(p->l, ray, nextFace, hit);
        }
        return flag;
    }
};

#endif  // !OBJECTKDTREE_H