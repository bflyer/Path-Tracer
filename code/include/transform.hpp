#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <vecmath.h>
#include "object3d.hpp"

// Ref: ver.2020

// transforms a 3D point using a matrix, returning a 3D point(对坐标点，添加平移效果)
static Vector3f transformPoint(const Matrix4f &mat, const Vector3f &point) {
    return (mat * Vector4f(point, 1)).xyz();
}

// transform a 3D direction using a matrix, returning a direction（对方向，不添加平移效果）
static Vector3f transformDirection(const Matrix4f &mat, const Vector3f &dir) {
    return (mat * Vector4f(dir, 0)).xyz();
}

class Transform : public Object3D {
public:
    Transform() {}

    // TODO: 优化掉 Object3D(obj->getMaterial()
    Transform(const Matrix4f &m, Object3D *obj) : o(obj), area(1), Object3D(obj->getMaterial()) {
        transform = m.inverse();
        bounds[0] = transformPoint(m, o->min());
        bounds[1] = transformPoint(m, o->max());
        bounds[2] = transformPoint(m, o->center());
    }

    ~Transform() {
    }

    virtual bool intersect(const Ray &r, Hit &h, float tmin) {
        // 1. 将仿射变换逆变换作用到射线上，得到一根新的光线
        Vector3f trSource = transformPoint(transform, r.getOrigin());
        Vector3f trDirection = transformDirection(transform, r.getDirection());
        Ray tr(trSource, trDirection);
        // 2. 用新的光线和物体求交
        bool inter = o->intersect(tr, h, tmin);
        if (inter) {
            // 3. 由 (M^(-1))^T 得到仿射变换后的法向
            h.set(h.getT(), h.getMaterial(), 
            transformDirection(transform.transposed(), h.getNormal()).normalized(),
            h.getColor(), h.getT() * r.getDirection() + r.getOrigin());
        }
        return inter;
    }

    Vector3f min() const override { return bounds[0]; }
    Vector3f max() const override { return bounds[1]; }
    Vector3f center() const override { return bounds[2]; }
    vector<Object3D *> getFaces() override { return {(Object3D *)this}; }
    double getArea() const override { return area; }

protected:
    Object3D *o; // un-transformed object
    Matrix4f transform;
    Vector3f bounds[3];
    double area;
};

#endif //TRANSFORM_H
