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
#endif // !UTILS_H