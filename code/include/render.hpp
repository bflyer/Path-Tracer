#ifndef PATH_TRACER_H
#define PATH_TRACER_H

#include <cmath>
#include <cstring>
#include <iostream>
#include <string>
#include <time.h>

#include "camera.hpp"
#include "constants.h"
#include "group.hpp"
#include "hit.hpp"
#include "hit_kdtree.hpp"
#include "image.hpp"
#include "light.hpp"
#include "ray.hpp"
#include "scene_parser.hpp"
#include "utils.hpp"
using namespace std;

// Whitted Style Path Tracing
// Ref: ver.2020
static Vector3f wtColor(Ray ray, const SceneParser& sceneParser, int depth = 0, int E = 1) {
    Group* group = sceneParser.getGroup();
    Vector3f color = Vector3f::ZERO;
    Vector3f opacity = Vector3f(1, 1, 1);

    // 如果深度超过阈值或者不透明度过小，则返回
    if (++depth > TRACE_DEPTH || opacity.max() < OPACITY_THRESHOLD) 
        return color;

    Hit hit;
    // 如果没有交点，直接返回背景色
    if (!group->intersect(ray, hit, TMIN)) {
        color += sceneParser.getBackgroundColor();
        return color;
    }

    // 交点属性
    Vector3f hitPos(ray.getOrigin() + ray.getDirection() * hit.getT());
    Material* material = hit.getMaterial();
    Vector3f N(hit.getNormal());
    Vector3f unitDir(ray.getDirection().normalized());

    // 1. 将颜色初始化为 ambient color（环境光），此处只考虑是否为光源
    color += material->getEmission();

    // 2. 计算每个光源对该点的贡献
    for (int li = 0; li < sceneParser.getNumLights(); li++) {
        Light* light = sceneParser.getLight(li);
        Vector3f L, lightColor;
        // 获得光照强度
        light->getIllumination(ray.pointAtParameter(hit.getT()), L, lightColor);
        // 计算局部光强
        color += hit.getMaterial()->Shade(ray, hit, L, lightColor);
    }

    // TODO: 修改为 hit.getColor()，并修改物体类
    // 不透明度衰减
    opacity = opacity * material->getDiffuseColor();
    
    // 3.A 若该表面是反射面(Reflection)
    if (material->getType().y() == 1) {
        // TODO: r2 =?= r1 + 2 * N 
        float cosine = Vector3f::dot(unitDir, N);
        Vector3f newDir = (ray.direction - N * (cosine * 2)).normalized();
        // color += material->getSpecular() * wtColor(Ray(hitPos, newDir), sceneParser, depth);
        color += wtColor(Ray(hitPos, newDir), sceneParser, depth);
    } 
    // 3.B 若该表面是透射面(Transmission)
    else if (material->getType().z() > 0) {
        // 折射率
        float n = material->getRefractRate();
        // 入射角余弦值
        float cos1 = Vector3f::dot(unitDir, N);
        // 如果光线来自物体内部，将法线反向
        if (cos1 > 0) {
            N.negate();
        }
        // 相对折射率
        n = 1 / n;    
        // 折射角余弦值
        float cos2 = sqrt(1 - n * n * (1 - cos1 * cos1));
        float R0 = ((1.0 - n) * (1.0 - n)) / ((1.0 + n) * (1.0 + n));
        float Rprob = R0 + (1.0 - R0) * pow((1.0 - cos1), 5.0);  // 反射概率（Schlick-approximation）
        
        Vector3f newDir;
        // 非全反射，粗略地只考虑折射
        // TODO：改成 R.R.
        if (cos2 > 0) {
            newDir = ((-n) * unitDir + (n * Vector3f::dot(unitDir, N) - cos2) * N).normalized();
            color += wtColor(Ray(hitPos, newDir), sceneParser, depth);
        }
        // 全反射
        else {
            float cosine = Vector3f::dot(ray.direction, N);
            newDir = (ray.direction - N * (cosine * 2)).normalized();
            // color += material->getSpecular() * wtColor(Ray(hitPos, newDir), sceneParser, depth);
            color += wtColor(Ray(hitPos, newDir), sceneParser, depth);
        }
    }
    // 4. 返回颜色
    return color;
}

// // Whitted Style Path Tracing-2020
// static Vector3f wtColor(Ray ray, const SceneParser& sceneParser, int depth = 0, int E = 0) {
//     Group* group = sceneParser.getGroup();
//     // cf 代表光线的贡献因子（颜色滤镜），随着光线的反射或折射逐渐衰减
//     Vector3f color(0, 0, 0), cf(1, 1, 1);
//     while (true) {
//         if (++depth > TRACE_DEPTH || cf.max() < 1e-3) return color;
//         // 判断camRay是否和场景有交点,返回最近交点的数据,存储在hit中.
//         Hit hit;
//         if (!group->intersect(ray, hit, TMIN)) {
//             color += sceneParser.getBackgroundColor();
//             return color;
//         }

//         // Path Tracing
//         // 获取下一个交点及相关性质
//         ray.origin += ray.direction * hit.getT();
//         Material* material = hit.getMaterial();
//         Vector3f refColor(material->getDiffuseColor()), N(hit.getNormal());

//         // Emission
//         color += material->getEmission() * cf;
//         cf = cf * refColor;
//         float type = RAND2;
//         if (type <= material->getType().x()) {  // diffuse
//             ray.direction = diffDir(N);
//         } else if (type <= material->getType().x() + material->getType().y()) {  // specular
//             float cost = Vector3f::dot(ray.direction, N);
//             ray.direction = (ray.direction - N * (cost * 2)).normalized();
//         } else {  // refraction
//             float n = material->getRefractRate();
//             float R0 = ((1.0 - n) * (1.0 - n)) / ((1.0 + n) * (1.0 + n));
//             if (Vector3f::dot(N, ray.direction) > 0) {  // inside the medium
//                 N.negate();
//                 n = 1 / n;
//             }
//             n = 1 / n;
//             float cost1 = -Vector3f::dot(N, ray.direction);  // 入射角 cosine
//             float cost2 =
//                 1.0 - n * n * (1.0 - cost1 * cost1);  // 折射角 cosine
//             float Rprob = R0 + (1.0 - R0) * pow(1.0 - cost1,
//                                                 5.0);  // Schlick-approximation
//             if (cost2 > 0 && RAND2 > Rprob) {           // refraction direction
//                 ray.direction =
//                     ((ray.direction * n) + (N * (n * cost1 - sqrt(cost2))))
//                         .normalized();
//             } else {  // reflection direction
//                 ray.direction = (ray.direction + N * (cost1 * 2));
//             }
//         }
//     }
// }

// Path Tracing
// Ref: smallpt
static Vector3f ptColor(Ray ray, const SceneParser& sceneParser, int depth = 0, int E = 1) {
    Group* group = sceneParser.getGroup();
    Vector3f color = Vector3f::ZERO;
    Hit hit;
    // 1. 求交：如果没有交点，直接返回背景色
    // TODO: 返回背景色还是返回黑色呢？
    if (!group->intersect(ray, hit, TMIN)) {
        return sceneParser.getBackgroundColor();
    }

    // 交点坐标
    Vector3f hitPos(ray.getOrigin() + ray.getDirection() * hit.getT());
    Material* material = hit.getMaterial();          // the hit object
    Vector3f f(material->getDiffuseColor());         // BRDF
    Vector3f n(hit.getNormal());
    Vector3f nl = Vector3f::dot(n, ray.getDirection()) < 0 ? n : -n;  // ensure the normal points outward

    // TODO：丢弃时选择 0 还是选择材质的 emission
    // 2. R.R.(Russian Roulette)
    if (++depth > RR_DEPTH) {
        if (RAND2 < RR_PROBABILITY)
            f = f * (1.0 / RR_PROBABILITY);
        else 
            return material->getEmission() * E;
    }

    float type = RAND2;
    // 3. 单独处理理想漫反射和理想镜面反射
    // Ideal DIFFUSE reflection(理想漫反射)
    if (material->getType().x() == 1){             
        // 记极角(polar angle) 为 theta，记方位角(Azimuthal angle) 为 phi    
        double phi = 2 * M_PI * RAND2;  
        double sinTheta2 = RAND2;              // sin(theta) ^ 2
        double sinTheta = sqrt(sinTheta2);     // sin(theta)
        
        // 构建局部坐标系
        Vector3f w = nl;  // 表面的单位法线
        // 通过判断法线向量 w 的 x 分量绝对值是否大于 0.1 来决定
        // 构造正交基的方式，以避免因法线几乎平行于 x 轴而导致的
        // 除以零问题。
        Vector3f u = (Vector3f::cross((fabs(w.x()) > .1 ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0)), w)).normalized();
        Vector3f v = Vector3f::cross(w, u);
        // 利用随机数 phi 和 sinTheta，以及正交基 u 和 v，
        // 通过球坐标系转换为笛卡尔坐标系的方式计算出新的反射方向向量 d
        Vector3f d = (u * cos(phi) * sinTheta + v * sin(phi) * sinTheta + w * sqrt(1 - sinTheta2)).normalized();
        float cosHit = Vector3f::dot(d, n);
        float c = (cosHit > 0 ? cosHit : -cosHit) * 2 * M_PI;
        // float c = (cosHit > 0 ? cosHit : -cosHit);

        // 对光源采样
        // Loop over any lights
        Vector3f e;
        Hit h;
        // for (Sphere* eObj : group->getEmissionObjList()){
        //     // 用 Realistic Ray Tracing 创建打向球体的随机光线方向
        //     Vector3f sw = eObj->getCenter() - hitPos;         // 发光球体球心与交点的连线
        //     Vector3f su = Vector3f::cross((fabs(sw.x()) > .1 ? Vector3f(0, 1, 1) : Vector3f(1, 0, 0)), sw).normalized();
        //     Vector3f sv = Vector3f::cross(sw, su);

        //     double cos_a_max = sqrt(1 - eObj->getRadius() * eObj->getRadius() / Vector3f::dot(hitPos - eObj->getCenter(), hitPos - eObj->getCenter()));
        //     double eps = RAND2;
        //     double cos_a = 1 - eps + eps * cos_a_max;
        //     double sin_a = sqrt(1 - cos_a * cos_a);
        //     double phi = 2 * M_PI * RAND2;
        //     Vector3f l = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;

        //     // Shoot shadow ray
        //     // if (group->intersect(Ray(hitPos, l.normalized()), h, TMIN) && h.getT() == l.length()){  // Check for occlusion with shadow ray
        //     if (group->intersect(Ray(hitPos, l.normalized()), h, TMIN)){  // shadow ray
        //         double omega = 2 * M_PI * (1 - cos_a_max);
        //         e = e + f * (eObj->getMaterial()->getEmission() * Vector3f::dot(l, nl) * omega) * M_1_PI;  // 1/pi for brdf
        //     }
        // }

        // for (Sphere* eObj : group->getEmissionObjList()){
        //     // 用 Realistic Ray Tracing 创建打向球体的随机光线方向
        //     Vector3f sw = eObj->getCenter() - hitPos;         // 发光球体球心与交点的连线
        //     Vector3f su = Vector3f::cross((fabs(sw.x()) > .1 ? Vector3f(0, 1, 1) : Vector3f(1, 0, 0)), sw).normalized();
        //     Vector3f sv = Vector3f::cross(sw, su);

        //     double cos_a_max = sqrt(1 - eObj->getRadius() * eObj->getRadius() / Vector3f::dot(hitPos - eObj->getCenter(), hitPos - eObj->getCenter()));
        //     double eps = RAND2;
        //     double cos_a = eps * cos_a_max;           // 在可视锥体内选一个极角
        //     double sin_a = sqrt(1 - cos_a * cos_a);   
        //     double phi = 2 * M_PI * RAND2;
        //     Vector3f l = (su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a).normalized();

        //     // Shoot shadow ray
        //     // if (group->intersect(Ray(hitPos, l.normalized()), h, TMIN) && h.getT() == l.length()){  // Check for occlusion with shadow ray
        //     if (group->intersect(Ray(hitPos, -l), h, TMIN)){  // shadow ray
        //         float cos1 = Vector3f::dot(-l, n);          // 光线与原交点法线的夹角余弦
        //         cos1 = cos1 > 0 ? cos1 : -cos1;
        //         float cos2 = Vector3f::dot(h.getNormal(), l);
        //         float invD2 = 1 / (h.getT() * h.getT());
        //         double sin_a_max = sqrt(1 - cos_a_max * cos_a_max);   // 最大半角正弦
        //         // float invArea = 1 / (2 * M_PI * eObj->getRadius() * (1 - sin_a_max));
        //         float invArea = 1 / (eObj->getRadius() * (1 - sin_a_max));
        //         e = e + f * eObj->getMaterial()->getEmission() * cos2 * invD2 * invArea;
        //         // e = e + f * eObj->getMaterial()->getEmission();
        //     }
        // }

        // return material->getEmission() * E + e + f * c * (ptColor(Ray(hitPos, d), sceneParser, depth, 1));
        return material->getEmission() + f *(ptColor(Ray(hitPos, d), sceneParser, depth, 1));
    }

    // Ideal SPECULAR reflection(理想镜面反射)
    else if (material->getType().y() == 1) {
        Vector3f d = ray.getDirection() - n * 2 * Vector3f::dot(ray.getDirection(), n);
        return material->getEmission() + f * (ptColor(Ray(hitPos, d), sceneParser, depth));
    }
        
    // Ideal dielectric REFRACTION(理想介质折射)
    // 反射光线初始化，直接套用反射公式
    Vector3f d = ray.getDirection() - n * 2 * Vector3f::dot(ray.getDirection(), n);
    Ray reflRay(hitPos, d); 
    // into = true: 光线从外而内；into = false：光线从内而外    
    bool into = Vector3f::dot(n, nl) > 0;                
    // nc: 外部介质折射率（如空气）
    // nt: 内部介质的折射率（如玻璃）
    double nc = 1, nt = material->getRefractRate();
    // nnt: 根据光线进出方向计算相对折射率
    double nnt = into ? nc/nt : nt/nc;
    // ddn: 入射角余弦值
    double ddn = Vector3f::dot(ray.getDirection(), nl);
    double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
    
    // Total internal reflection(全内反射检查)
    if (cos2t < 0) {
        return material->getEmission() + f *(ptColor(reflRay, sceneParser, depth));
    }
        
    // 计算折射方向
    Vector3f tdir = (ray.getDirection() * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalized();

    // 计算菲涅尔反射和折射系数
    // R0: 菲涅耳反射系数，最终反射率
    // Re: 最终折射率；Tr: 最终反射率
    // c: 辅助变量，用于确定光线能量在反射和折射间的分配
    double a = nt - nc, b = nt + nc, R0 = a * a / (b*b), c = 1 - (into ? -ddn : Vector3f::dot(tdir, n));
    double Re = R0 + (1 - R0) *c*c*c*c*c;
    double Tr = 1 - Re;
    // P: 用于决定追踪反射光线还是折射光线
    // RP: 追踪反射光线的权重；TP: 追踪折射光线的权重
    double P = .25 + .5*Re, RP = Re/P, TP = Tr / (1-P);

    // 递归深度大于阈值时使用俄罗斯轮盘赌
    return material->getEmission() + f * (depth > 2 ? (RAND2 < P ?
        ptColor(reflRay, sceneParser, depth) * RP : ptColor(Ray(hitPos, tdir), sceneParser, depth) * TP) :
        ptColor(reflRay, sceneParser, depth) * Re + ptColor(Ray(hitPos, tdir), sceneParser, depth) * Tr);
}

// // Path Tracing（递推版）
// // Ref: smallpt
// static Vector3f ptColor(Ray ray, const SceneParser& sceneParser, int depth = 0, int E = 1) {
//     Group* group = sceneParser.getGroup();
//     Vector3f color = Vector3f::ZERO;
//     Vector3f opacity = Vector3f(1, 1, 1);
//     Hit hit;
//     Ray r = ray;

//     while (true) {
//         // 1. 求交：如果没有交点，直接返回背景色
//         // TODO: 返回背景色还是返回黑色呢？
//         if (!group->intersect(r, hit, TMIN)) {
//             return color;
//         }
        
//         // 交点坐标
//         Vector3f hitPos(r.getOrigin() + r.getDirection() * hit.getT());
//         Material* material = hit.getMaterial();          // the hit object
//         Vector3f f(material->getDiffuseColor());         // BRDF
//         Vector3f n(hit.getNormal());
//         Vector3f nl = Vector3f::dot(n, r.getDirection()) < 0 ? n : -n;  // ensure the normal points outward
        
//         // TODO：丢弃时选择 0 还是选择材质的 emission
//         // 2. R.R.(Russian Roulette)
//         if (++depth > RR_DEPTH) {
//             if (RAND2 < RR_PROBABILITY)
//                 f = f * (1.0 / RR_PROBABILITY);
//             else 
//                 return color;
//         }

//         color = color + opacity * material->getEmission();
//         opacity = opacity * f;

//         float type = RAND2;
//         // 3. 单独处理理想漫反射和理想镜面反射
//         // Ideal DIFFUSE reflection(理想漫反射)
//         if (material->getType().x() == 1){             
//             // 记极角(polar angle) 为 theta，记方位角(Azimuthal angle) 为 phi    
//             double phi = 2 * M_PI * RAND2;  
//             double sinTheta2 = RAND2;              // sin(theta) ^ 2
//             double sinTheta = sqrt(sinTheta2);     // sin(theta)
            
//             // 构建局部坐标系
//             Vector3f w = nl;  // 表面的单位法线
//             // 通过判断法线向量 w 的 x 分量绝对值是否大于 0.1 来决定
//             // 构造正交基的方式，以避免因法线几乎平行于 x 轴而导致的
//             // 除以零问题。
//             Vector3f u = (Vector3f::cross((fabs(w.x()) > .1 ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0)), w)).normalized();
//             Vector3f v = Vector3f::cross(w, u);
//             // 利用随机数 phi 和 sinTheta，以及正交基 u 和 v，
//             // 通过球坐标系转换为笛卡尔坐标系的方式计算出新的反射方向向量 d
//             Vector3f d = (u * cos(phi) * sinTheta + v * sin(phi) * sinTheta + w * sqrt(1 - sinTheta2)).normalized();

//             r = Ray(hitPos, d);
//             continue;
//         }

//         // Ideal SPECULAR reflection(理想镜面反射)
//         else if (material->getType().y() == 1) {
//             r = Ray(hitPos, r.getDirection() - n * 2 * Vector3f::dot(r.getDirection(), n));
//             continue;
//         }
            
//         // Ideal dielectric REFRACTION(理想介质折射)
//         // 反射光线初始化，直接套用反射公式
//         Ray reflRay(hitPos, ray.getDirection() - n * 2 * Vector3f::dot(ray.getDirection(), n)); 
//         // into = true: 光线从外而内；into = false：光线从内而外    
//         bool into = Vector3f::dot(n, nl) > 0;                
//         // nc: 外部介质折射率（如空气）
//         // nt: 内部介质的折射率（如玻璃）
//         double nc = 1, nt = material->getRefractRate();
//         // nnt: 根据光线进出方向计算相对折射率
//         double nnt = into ? nc/nt : nt/nc;
//         // ddn: 入射角余弦值
//         double ddn = Vector3f::dot(ray.getDirection(), nl);
//         double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
//         // Total internal reflection(全内反射检查)
//         if (cos2t < 0) {
//             r = reflRay;
//             continue;
//         }
            
//         // 计算折射方向
//         Vector3f tdir = (ray.getDirection() * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalized();
//         // 计算菲涅尔反射和折射系数
//         // R0: 菲涅耳反射系数，最终反射率
//         // Re: 最终折射率；Tr: 最终反射率
//         // c: 辅助变量，用于确定光线能量在反射和折射间的分配
//         double a = nt - nc, b = nt + nc, R0 = a * a / (b*b), c = 1 - (into ? -ddn : Vector3f::dot(tdir, n));
//         double Re = R0 + (1 - R0) *c*c*c*c*c;
//         double Tr = 1 - Re;
//         // P: 用于决定追踪反射光线还是折射光线
//         // RP: 追踪反射光线的权重；TP: 追踪折射光线的权重
//         double P = .25 + .5*Re, RP = Re/P, TP = Tr / (1-P);
//         // R.R.
//         // 递归深度大于阈值时使用俄罗斯轮盘赌
//         if (RAND2 < P) {
//             opacity = opacity * RP;
//             r = reflRay;
//         } else {
//             opacity = opacity * TP;
//             r = Ray(hitPos, tdir);
//         }
//     }
// }


// Ray Casting
// Ref: ver.2020
static Vector3f rcColor(Ray ray, const SceneParser& sceneParser, int depth = 0, int E = 1) {
    Group* group = sceneParser.getGroup();
    Hit hit;
    Vector3f color = Vector3f::ZERO;
    
    // 如果没有交点，直接返回背景色
    // TODO: 确定阈值
    if (!group->intersect(ray, hit, 0)) {
        color += sceneParser.getBackgroundColor();
        return color;
    }
    // 如果有交点，计算每个光源提供的贡献
    for (int li = 0; li < sceneParser.getNumLights(); li++) {
        Light* light = sceneParser.getLight(li);
        Vector3f L, lightColor;
        // 获得光照强度
        light->getIllumination(ray.pointAtParameter(hit.getT()), L, lightColor);
        // 计算局部光强
        color += hit.getMaterial()->Shade(ray, hit, L, lightColor);
    }
    return color;
}

// Ref: ver.2020
// Ref: smallpt
class PathTracer {
    public:
        const SceneParser& sceneParser;
        int samplesPerPixel;
        const char* fout;
        Vector3f (*radiance)(Ray ray, const SceneParser& sceneParser, int depth, int E);

        // Constructer
        PathTracer(const SceneParser& sceneParser, int samplesPerPixel, const char* fout, const char* method)
            : sceneParser(sceneParser), 
            samplesPerPixel(samplesPerPixel),
            fout(fout) {
                if (strcmp(method, "pt") == 0) 
                    radiance = ptColor;
                else if (strcmp(method, "rc") == 0)
                    radiance = rcColor;
                else if (strcmp(method, "wt") == 0)
                    radiance = wtColor;
                else {
                    cout << "Invalid method: " << method << endl;
                    exit(1);
                }
                cout << "Path Tracing Constructed Successfully." << endl;
            }

        void render() {
            Camera* camera = sceneParser.getCamera();
            int w = camera->getWidth(), h = camera->getHeight();
            cout << "Width: " << w << " Height: " << h << endl;
            Image outImg(w, h);
            time_t start = time(NULL);  

#pragma omp parallel for schedule(dynamic, 1)  // OpenMP
            // 超分辨率采样因子，例如 2x2，4x4 等；invss2 是其平方的倒数
            const int superSample = 2;
            const float invss2 = 1.0f / (superSample * superSample);
            for (int yy = 0; yy < h * superSample; ++yy) {
                for (int xx = 0; xx < w * superSample; ++xx) {
                    // 计算实际输出图像的像素位置
                    int x = xx / superSample;
                    int y = yy / superSample;
                    
                    Vector3f color = Vector3f::ZERO;
                    for (int s = 0; s < samplesPerPixel; ++s) {
                        // 使用抖动（jittering）在子像素区域内采样
                        float subpixelX = (xx % superSample + RAND2) / superSample;
                        float subpixelY = (yy % superSample + RAND2) / superSample;
                        Ray camRay = camera->generateRay(Vector2f(x + subpixelX, y + subpixelY));
                        color += radiance(camRay, sceneParser, 0, 1);
                    }
                    // 对超级采样区域内的颜色求平均
                    color *= invss2;
                    outImg.SetPixel(x, y, (color / samplesPerPixel) + outImg.GetPixel(x, y));
                }
            }
            outImg.SaveBMP(fout);
        }
};

#endif  // !PATH_TRACER_H