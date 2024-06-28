#ifndef MATERIAL_H
#define MATERIAL_H

#include <cassert>
#include <vecmath.h>

#include "ray.hpp"
#include "hit.hpp"
#include "texture.h"

#include <iostream>

#define clampCos(x) ((x > 1) ? 1 : ((x < -1) ? -1 : x))

// Ref: ver.2020
class Material {
public:
    explicit Material(const Vector3f &d_color, 
                      const Vector3f &s_color = Vector3f::ZERO, 
                      const Vector3f &e_color = Vector3f::ZERO, 
                      float s = 0, float r = 0, Vector3f t = Vector3f(1, 0, 0),
                      const char *texture_filename = "", const char *bump_filename = "",
                      bool g = false, float m = 0, float rough = 0, float ior = 0)
            : diffuseColor(d_color), 
            specularColor(s_color), 
            emission(e_color),
            shininess(s), 
            refractRate(r), 
            type(t),
            texture(texture_filename),
            bump(bump_filename),
            glossy(g),
            metallic(m),
            roughness(rough),
            ior(ior) {
                if (e_color != Vector3f::ZERO && d_color != Vector3f::ZERO) 
                    emission = e_color * d_color;
            }

    virtual ~Material() = default;

    virtual Vector3f getDiffuseColor() const {
        return diffuseColor;
    }

    Vector3f getColor(float u, float v) const {
        if (!texture.pic)
            return diffuseColor;
        else
            return texture.getColor(u, v);
    }

    Vector3f getSpecular() const { return specularColor; }
    Vector3f getEmission() const { return emission; }
    Vector3f getType() const { return type; }
    bool getGlossy() const { return glossy; }
    float getMetallic() const { return metallic; }
    float getRoughness() const { return roughness; }
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

    /*
        1. I is the incident view direction(入射光线的方向)
        
        2. N is the normal at the intersection point
        
        3. ior is the material refractive index
        
        4. kr is the amount of light reflected(计算后得到的反射光的比例)
    */
    float fresnel(const Vector3f &I, const Vector3f &N, float ior) {
            float kr = 0.0;
            float cosi = clampCos(Vector3f::dot(I, N));
            float etai = 1, etat = ior;
            if (cosi > 0) {  std::swap(etai, etat); }
            // Compute sini using Snell's law
            float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
            // Total internal reflection
            if (sint >= 1) {
                kr = 1;
            }
            else {
                float cost = sqrtf(std::max(0.f, 1 - sint * sint));
                cosi = fabsf(cosi);
                float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
                float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
                kr = (Rs * Rs + Rp * Rp) / 2;
            }
            return kr;
            // As a consequence of the conservation of energy, transmittance is given by:
            // kt = 1 - kr;
    }

    /*
        Glossy 材质
        1. wi: 入射方向
        2. w0: 出射方向
        3. N: 法线
    */
    Vector3f eval(const Vector3f &wi, const Vector3f &wo, const Vector3f &N){
        float cosalpha = Vector3f::dot(N, wo);
        if (cosalpha > 0.0f) {
            // calculate the contribution of Microfacet model
            float F, G, D;

            F = fresnel(wi, N, ior);

            auto G_function = [&](const float& roughness, const Vector3f& wi, const Vector3f& wo, const Vector3f& N)
            {
                float A_wi, A_wo;
                A_wi = (-1 + sqrt(1 + roughness * roughness * pow(tan(acos(Vector3f::dot(wi, N))), 2))) / 2;
                A_wo = (-1 + sqrt(1 + roughness * roughness * pow(tan(acos(Vector3f::dot(wo, N))), 2))) / 2;
                float divisor = (1 + A_wi + A_wo);
                if (divisor < 0.001)
                    return 1.f;
                else
                    return 1.0f / divisor;
            };
            G = G_function(roughness, -wi, wo, N);

            auto D_function = [&](const float& roughness, const Vector3f& h, const Vector3f& N)
            {
                float cos_sita = Vector3f::dot(h, N);
                float divisor = (M_PI * pow(1.0 + cos_sita * cos_sita * (roughness * roughness - 1), 2));
                if (divisor < 0.001)
                    return 1.f;
                else 
                    return (roughness * roughness) / divisor;
            };
            Vector3f h = (-wi + wo).normalized();
            D = D_function(roughness, h, N);

            // energy balance
            Vector3f diffuse = (Vector3f(1.0f) - F) * diffuseColor / M_PI;
            Vector3f specular;
            float divisor= ((4 * (Vector3f::dot(N, -wi)) * (Vector3f::dot(N, wo))));
            if (divisor < 0.001)
                specular= Vector3f(1);
            else
                specular = F *G * D / divisor;

            return diffuse+specular;
        }
        else
            return Vector3f(0.0f);
    }

protected:
    Vector3f diffuseColor;         // 漫反射系数（albedo）
    Vector3f specularColor;        // 镜面反射系数
    float shininess;               // 高光指数
    Vector3f emission;             // 自发光系数
    float refractRate;              // 折射率
    Vector3f type;                 // 材质类型(如漫反射、镜面反射、光泽等)
    Texture texture, bump;         // 纹理与凹凸贴图
    bool glossy;                   // 是否为 glossy 材质
    float metallic;               // 金属度 [0, 1]
    float roughness;              // 粗糙度 [0, 1]
    float ior;                    // 反射率
    float clamp(float x) { return std::max(float(0), x); }   // 截断函数
};


#endif // MATERIAL_H