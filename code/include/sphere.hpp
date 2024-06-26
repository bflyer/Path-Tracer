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
        : Object3D(material), ballCenter(center), radius(radius) {}

    ~Sphere() override = default;

    bool intersect(const Ray &r, Hit &h, float tmin) override {
        // 1. 计算由光源指向球心的向量 l
        Vector3f o(r.getOrigin()), dir(r.getDirection().normalized());
        Vector3f l(ballCenter - o);

        // 2. 计算【球心】到【光线所在直线】的投影点（垂足）到光线起点的距离 foot
        float foot = Vector3f::dot(l, dir);

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
        float d2 = l.squaredLength() - foot * foot;
        // 若光源在球体外部且和视线方向同侧，但球心距离视线超过半径，则不相交
        if (d2 > radius * radius && position) {
            return false;
        }

        // 5. 计算【投影点】到【光线与球面的交点】的距离 delta
        float delta = sqrt(radius * radius - d2);

        // 6. 求解光线与球面的交点
        float t = position ? foot - delta : foot + delta;
        if (t < tmin || t >= h.getT()) {
            return false;
        } else {
            // 获取交点坐标及法向量
            Vector3f point(o + dir * t);
            Vector3f normal((point - ballCenter).normalized());
            h.set(t, material, normal);
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

    // Ref: ver.2020
    Ray randomRay(int axis=-1, long long int seed=0) const override {
        float u = 2 * myRandom(axis, seed) - 1, v = 2 * myRandom(axis, seed) - 1;
        float r2 = u * u + v * v;
        while(r2 >= 1) {
            ++seed;
            u = 2 * myRandom(axis, seed) - 1;
            v = 2 * myRandom(axis, seed) - 1;
            r2 = u * u + v * v;
        }
        Vector3f dir(2*u*sqrtf(1-r2), 2*v*sqrt(1-r2),1-2*r2);
        dir.normalize();
        return Ray(ballCenter + radius * dir, dir);
    }

    Vector3f min() const override { return ballCenter - radius; }
    Vector3f max() const override { return ballCenter + radius; }
    Vector3f center() const override { return ballCenter; }
    std::vector<Object3D *> getFaces() override { return {(Object3D *)this}; } 

protected:
    Vector3f ballCenter;  // 球心
    float radius;     // 半径
};
#endif