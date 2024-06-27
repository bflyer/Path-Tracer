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

    // 不透明度衰减
    opacity = opacity * hit.getColor();
    
    // 3.A 若该表面是反射面(Reflection)
    if (material->getType().y() == 1) {
        // TODO: r2 =?= r1 + 2 * N 
        float cosine = Vector3f::dot(unitDir, N);
        Vector3f newDir = (ray.direction - N * (cosine * 2)).normalized();
        // color += material->getSpecular() * wtColor(Ray(hitPos, newDir), sceneParser, depth);
        color += wtColor(Ray(hitPos, newDir), sceneParser, depth);
    } 
    // 3.B 若该表面是透射面(Transmission)
    else if (material->getType().z() == 1) {
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

// // Path Tracing
// // Ref: smallpt
// static Vector3f ptColor(Ray ray, const SceneParser& sceneParser, int depth = 0, int E = 1) {
//     Group* group = sceneParser.getGroup();
//     Vector3f color = Vector3f::ZERO;
//     Hit hit;

//     // 1. 求交：如果没有交点，直接返回背景色
//     if (!group->intersect(ray, hit, TMIN)) {
//         return sceneParser.getBackgroundColor();
//     }

//     // 交点坐标
//     Vector3f hitPos(ray.getOrigin() + ray.getDirection() * hit.getT());
//     Material* material = hit.getMaterial();          // the hit object
//     Vector3f f(hit.getColor());         // BRDF
//     Vector3f n(hit.getNormal());
//     Vector3f nl = Vector3f::dot(n, ray.getDirection()) < 0 ? n : -n;  // ensure the normal points outward

//     float p = f.max();
//     // 2. R.R.(Russian Roulette)
//     if (++depth > RR_DEPTH || !p) {  // 大于阈值或者打到光源，开始 R.R.（达到光源一定会返回）
//         if (RAND2 < p)
//             f = f * (1.0 / p);
//         else {
//             return material->getEmission() * E;
//         }
//     }

//     // 3. 单独处理理想漫反射和理想镜面反射
//     // Ideal DIFFUSE reflection(理想漫反射)
//     if (material->getType().x() == 1){      
//         // 记极角(polar angle) 为 theta，记方位角(Azimuthal angle) 为 phi    
//         double phi = 2 * M_PI * RAND2;  
//         double sinTheta2 = RAND2;              // sin(theta) ^ 2
//         double sinTheta = sqrt(sinTheta2);     // sin(theta)
        
//         // 构建局部坐标系
//         Vector3f w = nl;  // 表面的单位法线
//         // 通过判断法线向量 w 的 x 分量绝对值是否大于 0.1 来决定
//         // 构造正交基的方式，以避免因法线几乎平行于 x 轴而导致的
//         // 除以零问题。
//         Vector3f u = (Vector3f::cross((fabs(w.x()) > .1 ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0)), w)).normalized();
//         Vector3f v = Vector3f::cross(w, u);
//         // 利用随机数 phi 和 sinTheta，以及正交基 u 和 v，
//         // 通过球坐标系转换为笛卡尔坐标系的方式计算出新的反射方向向量 d
//         Vector3f d = (u * cos(phi) * sinTheta + v * sin(phi) * sinTheta + w * sqrt(1 - sinTheta2)).normalized();
//         float cosHit = Vector3f::dot(d, n);
//         // float c = (cosHit > 0 ? cosHit : -cosHit) * 2 * M_PI;
//         float c = (cosHit > 0 ? cosHit : -cosHit);
//         // 对光源采样
//         // Loop over any lights
//         Vector3f e = Vector3f::ZERO;
//         Hit h1, h2;

//         // for (Sphere* eObj : group->getEmissionObjList()){
//         //     // 用 Realistic Ray Tracing 创建打向球体的随机光线方向
//         //     Vector3f sw = (eObj->getCenter() - hitPos).normalized();         // 交点指向发光球体球心的单位向量
//         //     Vector3f su = Vector3f::cross((fabs(sw.x()) > .1 ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0)), sw).normalized();
//         //     Vector3f sv = Vector3f::cross(sw, su);
//         //     double cos_a_max2 = 1 - eObj->getRadius() * eObj->getRadius() / Vector3f::dot(hitPos - eObj->getCenter(), hitPos - eObj->getCenter());
//         //     if (cos_a_max2 < 0) continue;
//         //     double cos_a_max = sqrt(cos_a_max2);
//         //     double eps = RAND2;
//         //     double cos_a = 1 - eps + eps * cos_a_max;  // 先用半角公式缩到半角，取随机，然后再倍乘回来
//         //     double sin_a = sqrt(1 - cos_a * cos_a);
//         //     // double sin_a_max = sqrt(eObj->getRadius() * eObj->getRadius() / Vector3f::dot(hitPos - eObj->getCenter(), hitPos - eObj->getCenter()));
//         //     // double cos_a_max = sqrt(1 - sin_a_max * sin_a_max);
//         //     // double sin_a = RAND2 * sin_a_max;
//         //     // double cos_a = sqrt(1 - sin_a * sin_a);            
//         //     double phi = 2 * M_PI * RAND2;
//         //     Vector3f l = (su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a).normalized();
//         //     // Shoot shadow ray
//         //     // if (group->intersect(Ray(hitPos, l.normalized()), h, TMIN){  // Check for occlusion with shadow ray
//         //     if (group->intersect(Ray(hitPos, l), h1, TMIN))
//         //         if (eObj->intersect(Ray(hitPos, l), h2, TMIN) && h1.getT() == h2.getT()){  // shadow ray
//         //             float cos2 = Vector3f::dot(-h2.getNormal(), sw.normalized());
//         //             double omega = 2 * M_PI * (1 - cos_a_max);
//         //             float cosine = Vector3f::dot(l, nl);
//         //             cosine = cosine > 0 ? cosine : 0;
//         //             e = e + f * (eObj->getMaterial()->getEmission() * cosine * omega) * cos2;  // 1/pi for brdf
//         //             // e = e + f * (eObj->getMaterial()->getEmission() * cosine * omega) * cos2 / 2 * M_1_PI;  // 1/pi for brdf
//         //         }
//         // }
//         // return material->getEmission() * E + e + f * c * (ptColor(Ray(hitPos, d), sceneParser, depth, 0));
//         return material->getEmission() + f * c * (ptColor(Ray(hitPos, d), sceneParser, depth, 1));
//     }
    
//     // Ideal SPECULAR reflection(理想镜面反射)
//     else if (material->getType().y() == 1) {
//         Vector3f d = ray.getDirection() - n * 2 * Vector3f::dot(ray.getDirection(), n);
//         return material->getEmission() + f * (ptColor(Ray(hitPos, d), sceneParser, depth));
//     }
        
//     // Ideal dielectric REFRACTION(理想介质折射)
//     // 反射光线初始化，直接套用反射公式
//     Vector3f d = ray.getDirection() - n * 2 * Vector3f::dot(ray.getDirection(), n);
//     Ray reflRay(hitPos, d); 
//     // into = true: 光线从外而内；into = false：光线从内而外    
//     bool into = Vector3f::dot(n, nl) > 0;                
//     // nc: 外部介质折射率（如空气）
//     // nt: 内部介质的折射率（如玻璃）
//     double nc = 1, nt = material->getRefractRate();
//     // nnt: 根据光线进出方向计算相对折射率
//     double nnt = into ? nc/nt : nt/nc;
//     // ddn: 入射角余弦值
//     double ddn = Vector3f::dot(ray.getDirection(), nl);
//     double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
    
//     // Total internal reflection(全内反射检查)
//     if (cos2t < 0) {
//         return material->getEmission() + f *(ptColor(reflRay, sceneParser, depth));
//     }
        
//     // 计算折射方向
//     Vector3f tdir = (ray.getDirection() * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalized();

//     // 计算菲涅尔反射和折射系数
//     // R0: 菲涅耳反射系数，最终反射率
//     // Re: 最终折射率；Tr: 最终反射率
//     // c: 辅助变量，用于确定光线能量在反射和折射间的分配
//     double a = nt - nc, b = nt + nc, R0 = a * a / (b*b), c = 1 - (into ? -ddn : Vector3f::dot(tdir, n));
//     double Re = R0 + (1 - R0) *c*c*c*c*c;
//     double Tr = 1 - Re;
//     // P: 用于决定追踪反射光线还是折射光线
//     // RP: 追踪反射光线的权重；TP: 追踪折射光线的权重
//     double P = .25 + .5*Re, RP = Re/P, TP = Tr / (1-P);

//     // 递归深度大于阈值时使用俄罗斯轮盘赌
//     return material->getEmission() + f * (depth > 2 ? (RAND2 < P ?
//         ptColor(reflRay, sceneParser, depth) * RP : ptColor(Ray(hitPos, tdir), sceneParser, depth) * TP) :
//         ptColor(reflRay, sceneParser, depth) * Re + ptColor(Ray(hitPos, tdir), sceneParser, depth) * Tr);
// }

static Vector3f ptColor(Ray ray, const SceneParser& sceneParser, int depth = 0, int E = 1) {
    Group* group = sceneParser.getGroup();
    Vector3f color = Vector3f::ZERO;
    Hit hit;

    // 1. 求交：如果没有交点，直接返回背景色
    if (!group->intersect(ray, hit, TMIN)) {
        return sceneParser.getBackgroundColor();
    }
    // 交点坐标
    Vector3f hitPos(ray.getOrigin() + ray.getDirection() * hit.getT());
    Material* material = hit.getMaterial();          // the hit object
    Vector3f f(hit.getColor());         // BRDF
    Vector3f n(hit.getNormal());
    Vector3f nl = Vector3f::dot(n, ray.getDirection()) < 0 ? n : -n;  // ensure the normal points outward

    float p = f.max();
    // 2. R.R.(Russian Roulette)
    if (++depth > RR_DEPTH || !p) {  // 大于阈值或者打到光源，开始 R.R.（达到光源一定会返回）
        if (RAND2 < p)
            f = f * (1.0 / p);
        else {
            return material->getEmission() * E;
        }
    }

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
        // float c = (cosHit > 0 ? cosHit : -cosHit);
        // 对光源采样
        // Loop over any lights
        Vector3f e = Vector3f::ZERO;
        Hit h;
        for (Object3D* eObj : group->getEmissionObjList()){
            double area = eObj->getArea();        // 1. 获取光源面积
            Vector3f samplePoint = eObj->sample();       // 2. 对光源采样
            Vector3f sampleLine = samplePoint - hitPos;  // 着色点采样点的连线
            Vector3f sampleDir = sampleLine.normalized();   // 3. 得到光线方向
            if (group->intersect(Ray(hitPos, sampleDir), h, TMIN) &&    // 4. 检测是否被遮挡
                abs(h.getT() - sampleLine.length()) < TMIN) {  
            // if (group->intersect(Ray(hitPos, sampleDir), h, TMIN)) {
                double cos1 = Vector3f::dot(sampleDir, nl);              // 5. 计算光线与着色点法向量余弦
                cos1 = cos1 > 0 ? cos1 : -cos1;                         // TODO: 这里是否需要取绝对值？
                double cos2 = Vector3f::dot(sampleDir, hit.getNormal());  // 6. 计算光线与交点法向量余弦
                cos2 = cos2 > 0 ? cos2 : -cos2;
                e += eObj->getMaterial()->getEmission() * f * area * cos1 * cos2 / sampleLine.squaredLength();  // 7. 计算光源的颜色
            }
        }        
        return material->getEmission() * E + e + f * (ptColor(Ray(hitPos, d), sceneParser, depth, 0));
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
                    // cout << "(" << xx << ", " << yy << ") ";
                    // 计算实际输出图像的像素位置
                    int x = xx / superSample;
                    int y = yy / superSample;
                    Vector3f color = Vector3f::ZERO;
                    for (int s = 0; s < samplesPerPixel; ++s) { 
                        // cout << "s" << s << " ";
                        // 使用抖动（jittering）在子像素区域内采样
                        float subpixelX = (xx % superSample + RAND2) / superSample;
                        float subpixelY = (yy % superSample + RAND2) / superSample;
                        Ray camRay = camera->generateRay(Vector2f(x + subpixelX, y + subpixelY));
                        // if (xx == 16 && yy == 153) {
                        //     color += radiance(camRay, sceneParser, 0, 2);
                        // }
                        color += radiance(camRay, sceneParser, 0, 1);
                    }
                    // 对超级采样区域内的颜色求平均
                    // color = clampVec(color / samplesPerPixel) * invss2;
                    // outImg.SetPixel(x, y, color + outImg.GetPixel(x, y));
                    color *= invss2;
                    outImg.SetPixel(x, y, (color / samplesPerPixel) + outImg.GetPixel(x, y));
                }
            }
            outImg.SaveBMP(fout);
        }
};

#endif  // !PATH_TRACER_H