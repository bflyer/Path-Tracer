#ifndef PLANE_H
#define PLANE_H

#include "object3d.hpp"
#include <vecmath.h>
#include <cmath>
#include <cfloat>

// function: ax+by+cz=d
// choose your representation , add more fields and fill in the functions
// Ref: ver.2023(借鉴了 uv 坐标计算以及 uaxis 轴处理的细节)
class Plane : public Object3D {
public:
    Plane(const Vector3f &normal, float d, Material *m)
        : Object3D(m), 
        d(d / normal.length()), 
        normal(normal.normalized()), 
        area(1){
            // 如果法向几乎和 y 轴平行，选 (1, 0, 0)
            if (normal.x() < 1e-3 && normal.z() < 1e-3) {
                this->main_tangent = Vector3f::RIGHT;
            } else {  // 否则在其他两个轴向，用 (0, 1, 0) 和法向叉乘
                this->main_tangent = Vector3f::cross(Vector3f::UP, normal);
                this->main_tangent.normalize();
            }
            uaxis = Vector3f::cross(main_tangent, normal).normalized();
        }

    ~Plane() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        Vector3f o = r.getOrigin();
        Vector3f dir = r.getDirection();
        float dir_len = dir.length();
        dir.normalize();
        float c = Vector3f::dot(dir, normal);
        // 待改
        if (fabs(c) < 1e-3)
            return false;
        float t = (d - Vector3f::dot(normal, o)) / Vector3f::dot(normal, dir);
        t = t / dir_len;
        if (t < tmin || t > h.getT()) {
            return false;
        }
        Vector3f next_o = o + r.getDirection() * t;
        float v = next_o.y();
        float u = Vector3f::dot(next_o - r.getDirection() * normal, main_tangent);
        Vector2f grad = Vector2f::ZERO;
        float f = material->getBump().getDisturb(u, v, grad);

        Vector3f new_normal = normal;

        if (!(f < 1e-4 && f > -1e-4)) {
            new_normal += main_tangent * grad[0];
            new_normal += uaxis * grad[1];
            new_normal.normalize();
        }

        float uu = 0, vv = 0;
        getUV(uu, vv, next_o);

        h.set(t, this->material, new_normal, material->getColor(uu, vv), next_o);
        return true;
    }

    void getUV(float &u, float &v, const Vector3f &p) {
        v = Vector3f::dot(p - d * normal, uaxis) / 100;
        u = Vector3f::dot(p - d * normal, main_tangent) / 100;
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
    double area;      // 面积
    Vector3f main_tangent;
};

#endif //PLANE_H
		



