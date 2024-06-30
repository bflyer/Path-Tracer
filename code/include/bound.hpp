#ifndef BOUND_H
#define BOUND_H
#include <vecmath.h>
#include <vector>
#include "constants.h"
#include "ray.hpp"
using std::vector;

// Original
class AABB {
public:
    AABB() {
        bounds[0] = Vector3f(INF);
        bounds[1] = Vector3f(-INF);
    } 

    AABB (const Vector3f& min, const Vector3f& max) {
        bounds[0] = min;
        bounds[1] = max;
    }
    
    void set (const Vector3f& min, const Vector3f& max) {
        bounds[0] = min;
        bounds[1] = max;
    }

    // // TODO：可改用 maxElements 和 minElements 代替
    // void updateBound(const Vector3f &vec) {
    //     for (int i = 0; i < 3; ++i) {
    //         bounds[0][i] = bounds[0][i] < vec[i] ? bounds[0][i] : vec[i];
    //         bounds[1][i] = bounds[1][i] < vec[i] ? vec[i] : bounds[1][i];
    //     }
    // }
    // TODO：可改用 maxElements 和 minElements 代替
    void updateBound(const Vector3f &vec) {
        bounds[0] = minE(bounds[0], vec);
        bounds[1] = maxE(bounds[1], vec);
    }
    
    // 还没进来就已出去则无交点
    bool intersect(const Ray &r, float &t_min) {
        Vector3f o(r.getOrigin()), invDir(1 / r.getDirection());
        Vector3f t0 = (bounds[0] - o) * invDir;
        Vector3f t1 = (bounds[1] - o) * invDir;
        double tmin = INF, tmax = -INF;
        if (t0.x() < t1.x()) {
            tmin = t0.x();
            tmax = t1.x();
        } else {
            tmin = t1.x();
            tmax = t0.x();
        }
        if (t0.y() < t1.y()) {
            tmin = tmin > t0.y() ? tmin : t0.y();
            tmax = tmax < t1.y() ? tmax : t1.y();
        } else {
            tmin = tmin > t1.y() ? tmin : t1.y();
            tmax = tmax < t0.y() ? tmax : t0.y();
        }
        if (t0.z() < t1.z()) {
            tmin = tmin > t0.z() ? tmin : t0.z();
            tmax = tmax < t1.z() ? tmax : t1.z();
        } else {
            tmin = tmin > t1.z() ? tmin : t1.z();
            tmax = tmax < t0.z() ? tmax : t0.z();
        }
        if (tmin > tmax || tmax < 0) return false;
        t_min = tmin;
        return true;
    }

    Vector3f bounds[2];  // min and max of x, y, z
};


#endif  // !BOUND_H