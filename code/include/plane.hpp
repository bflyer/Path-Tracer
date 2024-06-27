#ifndef PLANE_H
#define PLANE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>
#include <cfloat>

// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions

class Plane : public Object3D {
public:
    Plane() : normal(Vector3f::UP), d(0) {}

    Plane(const Vector3f &normal, float d, Material *m)
        : Object3D(m), d(d), normal(normal) {
            uaxis = Vector3f::cross(Vector3f::UP, normal);
        }

    ~Plane() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        // 1. 计算 t
        float t = (d - Vector3f::dot(normal, r.getOrigin())) / Vector3f::dot(normal, r.getDirection());
        // 2. 设置交点
        if (std::isnan(t)) {
            return false;
        }
        if (t < tmin || t > h.getT()) {
            return false;
        } else {
            auto normal_ = Vector3f::dot(normal, r.getOrigin()) - d > 0 ? normal : -normal;
            h.set(t, material, normal_.normalized());
            return true;
        }
    }

    // Ref: ver.2020
    void getUV(float &u, float &v, const Vector3f &p) {
        v = p.y();
        u = Vector3f::dot(p - d * normal, uaxis);
    }

    Vector3f getNormal(float u, float v) {
        Vector2f grad(0);
        float f = material->getBump().getDisturb(u, v, grad);
        if (fabs(f) < FLT_EPSILON) return normal;
        if (uaxis.squaredLength() < FLT_EPSILON) return normal;
        return Vector3f::cross(uaxis + normal * grad[0],
                               Vector3f::UP + normal * grad[1])
            .normalized();
    }
    Vector3f min() const override {
        return -INF * Vector3f(fabs(normal.x()) < 1 - FLT_EPSILON,
                               fabs(normal.y()) < 1 - FLT_EPSILON,
                               fabs(normal.z()) < 1 - FLT_EPSILON) +
               normal * d;
    }
    Vector3f max() const override {
        return INF * Vector3f(fabs(normal.x()) < 1 - FLT_EPSILON,
                              fabs(normal.y()) < 1 - FLT_EPSILON,
                              fabs(normal.z()) < 1 - FLT_EPSILON) +
               normal * d;
    }
    Vector3f center() const override { return normal * d; }
    std::vector<Object3D *> getFaces() override { return {(Object3D *)this}; }

protected:
    Vector3f normal, uaxis;  // 平面法向量
    float d;          // 隐式表示中的 d
};

#endif //PLANE_H
		

