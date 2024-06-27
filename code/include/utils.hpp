#ifndef UTILS_H
#define UTILS_H
#include <cmath>
#include <cstdlib>
#include <random>

#include <Vector3f.h>

// Ref: ver.2020
#define PI 3.1415926536
const int prime[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 
    113, 127, 131, 137, 139, 149, 151, 157, 
    163, 167, 173,179, 181, 191, 193, 197, 
    199, 211, 223, 227, 229};

// 随机数生成器
static std::mt19937 mersenneTwister;
static std::uniform_real_distribution<float> uniform;
#define RAND (2.0 * uniform(mersenneTwister) - 1.0)         // [-1, 1)
#define RAND2 (uniform(mersenneTwister))                    // [0, 1)

inline double clamp(double x) { 
    return x < 0 ? 0 : (x > 1 ? 1 : x); 
}

inline Vector3f clampVec(const Vector3f& v) {
    return Vector3f(clamp(v.x()), clamp(v.y()), clamp(v.z()));
}

// 浮点数转整数颜色值：伽马校正(Gamma correction)逆变换
// 缩放至 [0, 255] 的范围，四舍五入
inline int toInt(float x) { 
    return int(pow(clamp(x), 1 / 2.2) * 255 + .5); 
}

// Ref: ver.2020
// 构建正交基
inline void construct_orthonormal_basis (const Vector3f& v1, Vector3f& v2, Vector3f& v3) {
    // 当 v1 的 x > y 时，投影到 x-z 平面上
    if (std::abs(v1.x()) > std::abs(v1.y())) {
        float invLen = 1.f / sqrtf(v1.x() * v1.x() + v1.z() * v1.z());
        v2 = Vector3f(-v1.z() * invLen, 0.0f, v1.x() * invLen);
    } 
    // 否则，投影到 y-z 平面上
    else {
        float invLen = 1.0f / sqrtf(v1.y() * v1.y() + v1.z() * v1.z());
        v2 = Vector3f(0.0f, v1.z() * invLen, -v1.y() * invLen);
    }
    v3 = Vector3f::cross(v1, v2);
}

// 单位半球上均匀分布的随机向量
/*
    cosTheta: cosine of polar angle
    randPhi: random number for azimuthal angle
*/
inline Vector3f hemisphere(float cosTheta, float randPhi) {
    const float r = sqrt(1.0 - cosTheta * cosTheta);
    const float phi = 2 * PI * randPhi;

    const float x = r * cos(phi);
    const float y = r * sin(phi);

    return Vector3f(x, y, cosTheta);
}

// 遵循余弦权重分布（Lambert 余弦定律分布）的单位半球向量
inline Vector3f cosineHemisphere(float cosTheta, float randPhi) {
    const float r = sqrt(cosTheta);
    const float phi = 2 * PI * randPhi;

    const float x = r * cos(phi);
    const float y = r * sin(phi);

    return Vector3f(x, y, sqrt(std::max(0.0f, 1 - cosTheta)));
}

// 基于低分歧序列(Low Discrepancy Sequence)的伪随机数生成方法
// 有时候被用于准蒙特卡洛积分(QMC)等需要高质量随机样本的场景
inline float randomQMC(int axis, long long int seed) {
    int base = prime[axis];
    float f = 1, res = 0;
    while (seed > 0) {
        f /= base;
        res += f * (seed % base);
        seed /= base;
    }
    return res;
}

// 如果提供 axis 则使用 randomQMC，否则使用 RAND2
inline float myRandom(int axis=-1, long long int seed=0) {
    if (axis == -1) return RAND2;
    return randomQMC(axis, seed);
}

// 生成一个随机的漫反射方向
inline Vector3f diffDir (const Vector3f& normal, int depth=0, long long int seed=0) {
    Vector3f rotX, rotY;
    construct_orthonormal_basis(normal, rotX, rotY);
    return Matrix3f(rotX, rotY, normal) * cosineHemisphere(myRandom(2 * depth+1, seed), myRandom(2 * depth+2, seed));
}

#endif // !UTILS_H