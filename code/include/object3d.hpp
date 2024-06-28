#ifndef OBJECT3D_H
#define OBJECT3D_H

#include "ray.hpp"
#include "hit.hpp"
#include "material.hpp"
#include <vector>

// Base class for all 3d entities.
class Object3D {
public:
    Object3D() : material(nullptr), area(1) {}

    virtual ~Object3D() = default;

    explicit Object3D(Material *material) {
        this->material = material;
    }

    Material *getMaterial() const { return material; }

    // Intersect Ray with this object. If hit, store information in hit structure.
    virtual bool intersect(const Ray &r, Hit &h, float tmin) = 0;
    virtual std::vector<Object3D *> getFaces() { return {this}; }
    virtual Vector3f min() const { return Vector3f(); }
    virtual Vector3f max() const { return Vector3f(); }
    virtual Vector3f center() const { return Vector3f(); }
    virtual double getArea() const { return getArea(); }
    virtual Vector3f sample() const { return Vector3f(); }
    
protected:
    Material *material;
    double area;
};

#endif

