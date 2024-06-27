#ifndef SPHERE_H
#define SPHERE_H

#include "object3d.hpp"
#include "utils.hpp"
#include <vecmath.h>
#include <cmath>

class Sphere : public Object3D {
public:
    Sphere() : radius(0), ballCenter(Vector3f::ZERO) {}

    Sphere(const Vector3f &center, float radius, Material *material)
        : Object3D(material), ballCenter(center), radius(radius), area(4 * M_PI * radius * radius) {}

    ~Sphere() override = default;

    bool intersect(const Ray &r, Hit &h, double tmin) override {
        // 1. 计算由光源指向球心的向量 l
        Vector3f o(r.getOrigin()), dir(r.getDirection().normalized());
        Vector3f l(ballCenter - o);

        // 2. 计算【球心】到【光线所在直线】的投影点（垂足）到光线起点的距离 foot
        double foot = Vector3f::dot(l, dir);

        // 3. 判断光源是否位于球体内部
        bool position;  // false 表示合法内部/表面， true 表示合法外部
        if (l.length() < this->radius) {  // 内部
            position = false;
        } else if (l.length() == this->radius) {   // 表面
            if (foot <= 0) {              
                return false;
            } else {
                position = false;
            }
        } else {  // 外部
            // 若光源在球体外部并且和视线方向异侧，则必然不相交
            if (foot <= 0) {
                return false;
            } else {
                position = true;
            }
        }

        // 4. 计算【球心】到【光线所在直线】的距离的平方 d2
        double d2 = l.squaredLength() - foot * foot;
        // 若光源在球体外部且和视线方向同侧，但球心距离视线超过半径，则不相交
        if (d2 > radius * radius && position) {
            return false;
        }

        // 5. 计算【投影点】到【光线与球面的交点】的距离 delta
        double delta = sqrt(radius * radius - d2);

        // 6. 求解光线与球面的交点
        double t = position ? foot - delta : foot + delta;
        if (t < tmin || t >= h.getT()) {
            return false;
        } else {
            // 获取交点坐标及法向量
            Vector3f point(o + dir * t);
            Vector3f OP = point - ballCenter;
            Vector3f normal(OP.normalized());
            // Vector3f normal((ballCenter - point).normalized());
            float u = 0.5 + atan2(normal.x(), normal.z()) / (2 * M_PI);
            float v = 0.5 - asin(normal.y()) * M_1_PI;
            h.set(t, material, getNormal(normal, OP, u, v), 
                  material->getColor(u, v), o + dir * t);
            return true;
        }
    }

    const Vector3f& getCenter() const {
        return ballCenter;
    } 

    float getRadius() {
        return radius;
    }

    // Ref: ver.2020
    // 计算 u-v 坐标下的法向
    Vector3f getNormal(const Vector3f &n, const Vector3f &p, float u, float v) {
        Vector2f grad(0);
        float f = material->getBump().getDisturb(u, v, grad);
        if (fabs(f) < FLT_EPSILON) return n;
        float phi = u * 2 * M_PI, theta = M_PI - v * M_PI;
        Vector3f pu(-p.z(), 0, p.x()),
            pv(p.y() * cos(phi), -radius * sin(theta), p.y() * sin(phi));
        if (pu.squaredLength() < FLT_EPSILON) return n;
        return Vector3f::cross(pu + n * grad[0] / (2 * M_PI),
                               pv + n * grad[1] / M_PI)
            .normalized();
    }

    Vector3f sample() const override {
        Vector3f randomDir(RAND2, RAND2, RAND2);
        return ballCenter + randomDir.normalized() * radius;
    }

    Vector3f min() const override { return ballCenter - radius; }
    Vector3f max() const override { return ballCenter + radius; }
    Vector3f center() const override { return ballCenter; }
    std::vector<Object3D *> getFaces() override { return {(Object3D *)this}; } 
    double getArea() const override{ return area; }

protected:
    Vector3f ballCenter;  // 球心
    float radius;     // 半径
    double area;
};
#endif