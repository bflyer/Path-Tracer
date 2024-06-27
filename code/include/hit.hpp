#ifndef HIT_H
#define HIT_H

#include <vecmath.h>
#include "ray.hpp"
#include "constants.h"

// Ref: ver.2020
class Material;

class Hit {
public:
    // constructors
    Hit()
        : material(nullptr),
            t(INF),
            r2(INIT_RADIUS),
            attenuation(Vector3f(1.0)),
            normal(Vector3f::ZERO),
            flux(Vector3f::ZERO),
            fluxLight(Vector3f::ZERO),
            color(Vector3f::ZERO),
            dir(Vector3f::ZERO),
            p(Vector3f::ZERO) {}

    Hit(float _t, Material *m, const Vector3f &n)
        : material(m),
            t(_t),
            r2(INIT_RADIUS),
            attenuation(Vector3f(1.0)),
            normal(n),
            flux(Vector3f::ZERO),
            fluxLight(Vector3f::ZERO),
            color(Vector3f::ZERO),
            dir(Vector3f::ZERO),
            p(Vector3f::ZERO) {}

    Hit(const Hit &h)
        : material(h.material),
            t(h.t),
            r2(INIT_RADIUS),
            attenuation(Vector3f(1.0)),
            normal(h.normal),
            flux(Vector3f::ZERO),
            fluxLight(Vector3f::ZERO),
            color(Vector3f::ZERO),
            dir(Vector3f::ZERO),
            p(Vector3f::ZERO) {}

    // destructor
    ~Hit() = default;

    float getT() const {
        return t;
    }

    Material *getMaterial() const { return material; }

    const Vector3f &getNormal() const { return normal; }
    const Vector3f &getColor() const { return color; }

    void set(float _t, Material *m, const Vector3f &n, const Vector3f &c,
             const Vector3f &_p) {
        t = _t;
        material = m;
        normal = n;
        color = c;
        p = _p;
    }

private:
    float t, r2;                // 表示交点到光线起点的距离
    Material *material;     // 交点处的材质属性
    Vector3f normal;        // 交点处的法线
    Vector3f color, flux, fluxLight, attenuation;  // 颜色、光通量、光源发出的光通量，光强衰减因子
    Vector3f dir, p;        // 交点处的入射方向、交点位置
};

inline std::ostream &operator<<(std::ostream &os, const Hit &h) {
    os << "Hit <" << h.getT() << ", " << h.getNormal() << ">";
    return os;
}

#endif // HIT_H
