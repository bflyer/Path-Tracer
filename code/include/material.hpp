// #ifndef MATERIAL_H
// #define MATERIAL_H

// #include <cassert>
// #include <vecmath.h>

// #include "ray.hpp"
// #include "hit.hpp"
// #include "texture.h"

// #include <iostream>

// // Ref: ver.2020
// class Material {
// public:

//     explicit Material(const Vector3f &d_color, 
//                       const Vector3f &s_color = Vector3f::ZERO, 
//                       const Vector3f &e_color = Vector3f::ZERO, 
//                       float s = 0, float r = 0, Vector3f t = Vector3f(1, 0, 0),
//                       const char *texture_filename = "", const char *bump_filename = "")
//             : diffuseColor(d_color), 
//             specularColor(s_color), 
//             emission(e_color),
//             shininess(s), 
//             refractRate(r), 
//             type(t), 
//             texture(texture_filename), 
//             bump(bump_filename) {}

//     virtual ~Material() = default;

//     virtual Vector3f getDiffuseColor() const {
//         return diffuseColor;
//     }

//     Vector3f getColor(float u, float v) {
//         if (!texture.pic) 
//             return diffuseColor;
//         else
//             return texture.getColor(u, v);
//     }

//     Vector3f getEmission() const { return emission; }
//     Vector3f getType() const { return type; }
//     float getRefractRate() const { return refractRate; }

//     // 单个光源着色
//     Vector3f Shade(const Ray &ray, const Hit &hit,
//                    const Vector3f &dirToLight, const Vector3f &lightColor) {
//         Vector3f shaded = Vector3f::ZERO;

//         // Step 1: Compute Diffuse
//         Vector3f N = hit.getNormal();           // 交点处法向
//         Vector3f L = dirToLight.normalized();                // 像素到光源的方向
//         Vector3f diffuse = diffuseColor * clamp(Vector3f::dot(N, L));

//         // Step 2: Compute Specular
//         Vector3f V = (-1) * ray.getDirection().normalized();   // 交点到像素的方向
//         Vector3f R = (2 * Vector3f::dot(N, L) * N - L).normalized();   // 反射光
//         Vector3f specular = specularColor * std::pow(clamp(Vector3f::dot(V, R)), shininess);

//         // Step 3: Compute Ambient(Ignore for now)
//         // Step 4: Compute Total
//         shaded = lightColor * (diffuse + specular);
//         return shaded;
//     }

// protected:
//     Vector3f diffuseColor;         // 漫反射系数
//     Vector3f specularColor;        // 镜面反射系数
//     float shininess;               // 高光指数
//     Vector3f emission;             // 自发光系数
//     float refractRate;              // 折射率
//     Vector3f type;                 // 材质类型(如漫反射、镜面反射、光泽等)
//     Texture texture, bump;         // 纹理与凹凸贴图
//     float clamp(float x) { return std::max(float(0), x); }   // 截断函数
// };


// #endif // MATERIAL_H

#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>

#include "ray.hpp"
#include "hit.hpp"
#include "texture.h"

#include <iostream>

// Ref: ver.2020
class Material {
public:

    explicit Material(const Vector3f &d_color, 
                      const Vector3f &s_color = Vector3f::ZERO, 
                      const Vector3f &e_color = Vector3f::ZERO, 
                      float s = 0, float r = 0, Vector3f t = Vector3f(1, 0, 0),
                      const char *texture_filename = "", const char *bump_filename = "")
            : diffuseColor(d_color), 
            specularColor(s_color), 
            emission(e_color),
            shininess(s), 
            refractRate(r), 
            type(t),
            texture(texture_filename),
            bump(bump_filename) {}

    virtual ~Material() = default;

    virtual Vector3f getDiffuseColor() const {
        return diffuseColor;
    }

    Vector3f getColor(float u, float v) {
        return diffuseColor;
    }

    Vector3f getEmission() const { return emission; }
    Vector3f getType() const { return type; }
    float getRefractRate() const { return refractRate; }

    Texture getBump() const { return bump; }

    // 单个光源着色
    Vector3f Shade(const Ray &ray, const Hit &hit,
                   const Vector3f &dirToLight, const Vector3f &lightColor) {
        Vector3f shaded = Vector3f::ZERO;

        // Step 1: Compute Diffuse
        Vector3f N = hit.getNormal();           // 交点处法向
        Vector3f L = dirToLight.normalized();                // 像素到光源的方向
        Vector3f diffuse = diffuseColor * clamp(Vector3f::dot(N, L));

        // Step 2: Compute Specular
        Vector3f V = (-1) * ray.getDirection().normalized();   // 交点到像素的方向
        Vector3f R = (2 * Vector3f::dot(N, L) * N - L).normalized();   // 反射光
        Vector3f specular = specularColor * std::pow(clamp(Vector3f::dot(V, R)), shininess);

        // Step 3: Compute Ambient(Ignore for now)
        // Step 4: Compute Total
        shaded = lightColor * (diffuse + specular);
        return shaded;
    }

protected:
    Vector3f diffuseColor;         // 漫反射系数
    Vector3f specularColor;        // 镜面反射系数
    float shininess;               // 高光指数
    Vector3f emission;             // 自发光系数
    float refractRate;              // 折射率
    Vector3f type;                 // 材质类型(如漫反射、镜面反射、光泽等)
    Texture texture, bump;         // 纹理与凹凸贴图
    float clamp(float x) { return std::max(float(0), x); }   // 截断函数
};


#endif // MATERIAL_H