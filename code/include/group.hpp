#ifndef GROUP_H
#define GROUP_H

#include "object3d.hpp"
#include "sphere.hpp"
#include "ray.hpp"
#include "hit.hpp"
#include "obj_kdtree.hpp"

#include <iostream>
#include <vector>

class Group : public Object3D {

public:

    Group() {}

    explicit Group (const vector<Object3D *> &objs, const vector<Object3D *> &eObjs)
        : objList(objs),
        eObjList(eObjs){
            kdTree = new ObjectKDTree(&objList);
        } 

    ~Group() override {}

    // KD-Tree + 包围盒加速求最近交点
    bool intersect(const Ray &r, Hit &h, float tmin) override {
        return kdTree->intersect(r, h, TMIN);
    }

    // // 求最近交点
    // bool intersect(const Ray &r, Hit &h, double tmin) override {
    //     bool flag = false;
    //     for (auto obj : objList)
    //         if (obj) flag |= obj->intersect(r, h, tmin);
    //     return flag;
    // }

    // 求最近交点
    bool sequentialSearch(const Ray &r, Hit &h, float tmin) {
        bool flag = false;
        for (auto obj : objList)
            if (obj) flag |= obj->intersect(r, h, tmin);
        return flag;
    }
    
    int getGroupSize() {
        return objList.size();
    }

    const std::vector<Object3D*>& getEmissionObjList() const {
        return eObjList;
    }

private:
    ObjectKDTree *kdTree;
    std::vector<Object3D*> objList;
    std::vector<Object3D*> eObjList;
};

#endif