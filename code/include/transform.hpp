#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <vecmath.h>
#include "object3d.hpp"

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

    Transform(const Matrix4f &m, Object3D *obj) : o(obj) {
        transform = m.inverse();
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
            h.set(h.getT(), h.getMaterial(), transformDirection(transform.transposed(), h.getNormal()).normalized());
        }
        return inter;
    }

protected:
    Object3D *o; // un-transformed object
    Matrix4f transform;
};

#endif //TRANSFORM_H
