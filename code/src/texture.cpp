#define STB_IMAGE_IMPLEMENTATION
// Ref: ver.2020
#include "texture.h"

#include <iostream>

#include "stb_image.h"
using namespace std;
// TODO：为什么要加上 w 和 h 再取模？
// 获取索引
// int Texture::getIdx(float u, float v) const {
//     int pw = int(u * w + w) % w, ph = int(v * h + h) % h;
//     return ph * w * channel + pw * channel;
// }

// 获取灰度值
float Texture::getGray(int idx) const { return (pic[idx] / 255. - 0.5) * 2; }

Texture::Texture(const char *textureFile) {
    if (strlen(textureFile) > 0) {
        pic = stbi_load(textureFile, &w, &h, &channel, 0);
        printf("Texture file: %s loaded. Size: %dx%dx%d\n", textureFile, w, h,
               channel);
    } else {
        pic = nullptr;
    }
}

// Vector3f Texture::getColor(int idx) const {
//     // std::cout << pic[idx] << " " << pic[idx + 1] << " " << pic[idx + 2] << std::endl;
//     return Vector3f(pic[idx], pic[idx + 1], pic[idx + 2]) / 255.;
// }

// Vector3f Texture::getColor(int u, int v) const {
//     // 若图片为空，直接返回 0
//     if (!pic) return Vector3f::ZERO;
//     u = u > w - 1 ? w - 1 : u;
//     v = v > h - 1 ? h - 1 : v;
//     int idx = (v * w  + u) * channel;
//     // std::cout << "int: " << pic[idx] << " " << pic[idx + 1] << " " << pic[idx + 2] << std::endl;
//     // Vector3f temp= Vector3f(pic[idx], pic[idx + 1], pic[idx + 2]) / 255.;
//     // temp.print();
//     return Vector3f(pic[idx], pic[idx + 1], pic[idx + 2]) / 255.;
// }

// Vector3f Texture::getColor(float u, float v) const {
//     // 若图片为空则直接返回 0 
//     if (!pic) return Vector3f::ZERO;
    
//     // 规范到 [0,1] 区间并放大
//     u -= int(u);
//     v -= int(v);
//     u = u < 0 ? 1 + u : u;
//     v = v < 0 ? 1 + v : v;
//     u = u * w;
//     v = h * (1 - v);  // 注意 u-v 坐标的 v 从下往上，但数组索引从上往下
    
//     // 取整，获得当前纹理坐标最近的左下角像素的整数坐标
//     int iu = (int)u, iv = (int)v;
//     float alpha = u - iu, beta = v - iv;
    
//     // 双线性插值
//     Vector3f ret(0, 0, 0);
//     ret += (1 - alpha) * (1 - beta) * getColor(iu, iv);
//     ret += alpha * (1 - beta) * getColor(iu + 1, iv);
//     ret += (1 - alpha) * beta * getColor(iu, iv + 1);
//     ret += alpha * beta * getColor(iu + 1, iv + 1);
//     // if (ret.x() > 1 || ret.y() > 1 || ret.z() > 1 || ret.x() < 0 || ret.y() < 0 || ret.z() < 0)
//         // std::cout << "float: " << ret.x() << " " << ret.y() << " " << ret.z() << std::endl;
//     return ret;
//     // TODO：为直接索引设计分支，加速
//     // int idx = getIdx(u, v);
//     // return Vector3f(pic[idx], pic[idx + 1], pic[idx + 2]) / 255.;
// }

// 扰动纹理
float Texture::getDisturb(float u, float v, Vector2f &grad) const {
    if (!pic) return 0;
    float disturb = getGray(getIdx(u, v));
    float du = 1.0 / w, dv = 1.0 / h;
    // 分别计算水平和竖直方向梯度
    grad[0] = w * (getGray(getIdx(u + du, v)) - getGray(getIdx(u - du, v))) / 2.0;
    grad[1] = h * (getGray(getIdx(u, v + dv)) - getGray(getIdx(u, v - dv))) / 2.0;
    return disturb;
}

int Texture::getIdx(float u, float v) const {
     // 减去整数部分，得到小数部分
    u -= int(u);
    v -= int(v);
    // 将小数部分限制在 [0, 1) 范围内
    u = u < 0 ? 1 + u : u;
    v = v < 0 ? 1 + v : v;
    // 根据 u 和 v 计算在纹理图像中的像素坐标
    int x = u * w;
    int y = v * h;
    // 将像素坐标限制在图像边界内
    x = x < 0 ? 0 : x;
    x = x > w - 1 ? w - 1 : x;
    y = y < 0 ? 0 : y;
    y = y > h - 1 ? h - 1 : y;
    // 将像素坐标限制在图像边界内
    int index = (y * w + x) * channel;
    return index;
};

Vector3f Texture::getColor(float u, float v) const {
     // 将小数部分限制在 [0, 1) 范围内
    u -= int(u);
    v -= int(v);
    // 将小数部分限制在 [0, 1) 范围内
    u = u < 0 ? 1 + u : u;
    v = v < 0 ? 1 + v : v;
    // 使用双线性插值计算纹理颜色
    u = u * w;
    v = h * (1 - v);
    int iu = (int)u, iv = (int)v;
    Vector3f ret_color = Vector3f::ZERO;
    float s = u - iu;
    float t = v - iv;
    // 使用双线性插值计算纹理颜色
    Vector3f color1 = (1 - s) * getColor(iu, iv + 1) + s * getColor(iu + 1, iv + 1);
    Vector3f color2 = (1 - s) * getColor(iu, iv) + s * getColor(iu + 1, iv);
    ret_color += (1 - t) * color2;
    ret_color += t * color1;
    return ret_color;
}

Vector3f Texture::getColor(int x, int y) const {
    x = x < 0 ? 0 : x;
    x = x > w - 1 ? w - 1 : x;
    y = y < 0 ? 0 : y;
    y = y > h - 1 ? h - 1 : y;
    int index = (y * w + x) * channel;
    return Vector3f(pic[index + 0], pic[index + 1], pic[index + 2]) / 255.0;
}