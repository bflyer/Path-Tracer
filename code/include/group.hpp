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

    explicit Group (const vector<Object3D *> &objs, const vector<Sphere *> &eObjs)
        : objList(objs),
        eObjList(eObjs){
            kdTree = new ObjectKDTree(&objList);
        } 

    ~Group() override {}

    // // KD-Tree + 包围盒加速求最近交点
    // bool intersect(const Ray &r, Hit &h, double tmin) override {
    //     return kdTree->intersect(r, h, TMIN);
    // }

    // 求最近交点
    bool intersect(const Ray &r, Hit &h, double tmin) override {
        bool flag = false;
        for (auto obj : objList)
            if (obj) flag |= obj->intersect(r, h, tmin);
        return flag;
    }

    // 求最近交点
    bool sequentialSearch(const Ray &r, Hit &h, float tmin) {
        bool flag = false;
        for (auto obj : objList)
            if (obj) flag |= obj->intersect(r, h, tmin);
        return flag;
    }

    // TODO: 待删除
    // void addObject(int index, Object3D *obj) {
    //     // TODO: 为什么需要有索引呢？直接插在最后不行吗？
    //     objList.push_back(obj);
    //     // objList.insert(objList.begin() + index, obj);
    // }

    // void addEmissionObject(Sphere *obj) {
    //     emissionObjList.push_back(obj);
    // }
    
    int getGroupSize() {
        return objList.size();
    }

    const std::vector<Sphere*>& getEmissionObjList() const {
        return eObjList;
    }

private:
    ObjectKDTree *kdTree;
    std::vector<Object3D*> objList;
    std::vector<Sphere*> eObjList;
};

#endif