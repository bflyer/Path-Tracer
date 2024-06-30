#define STB_IMAGE_IMPLEMENTATION
#define NORMALIZE(x) ((x < 0) ? 1 + x : x)  // 将小于 0 的小数置于 [0, 1] 范围内
#define NORMALIZE2(x, y) ((x < 0) ? 0 : ((x > y - 1) ? (y - 1) : x))  // 将数值限制在 [0, y - 1] 内

// Ref: ver.2023
#include "texture.h"
#include "stb_image.h"

#include <iostream>
using namespace std;

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

// 添加扰动
float Texture::getDisturb(float u, float v, Vector2f &grad) const {
    if (!pic) return 0;
    float disturb = getGray(getIdx(u, v));
    float du = 1.0 / w, dv = 1.0 / h;
    // 分别计算水平和竖直方向梯度
    grad[0] = w * (getGray(getIdx(u + du, v)) - getGray(getIdx(u - du, v))) / 2.0;
    grad[1] = h * (getGray(getIdx(u, v + dv)) - getGray(getIdx(u, v - dv))) / 2.0;
    return disturb;
}

// 获取在图片中对应的坐标索引
int Texture::getIdx(float u, float v) const {
     // 减去整数部分，得到小数部分
    u -= int(u);
    v -= int(v);
    // 将小数部分限制在 [0, 1) 范围内
    u = NORMALIZE(u);
    v = NORMALIZE(v);
    // 根据 u 和 v 计算在纹理图像中的像素坐标
    int x = u * w;
    int y = v * h;
    // 将像素坐标限制在图像边界内
    x = NORMALIZE2(x, w);
    y = NORMALIZE2(y, h);
    // 将像素坐标限制在图像边界内
    int index = (y * w + x) * channel;
    return index;
};

// 获取浮点坐标点的颜色
Vector3f Texture::getColor(float u, float v) const {
     // 将小数部分限制在 [0, 1) 范围内
    u -= int(u);
    v -= int(v);

    // 将小数部分限制在 [0, 1) 范围内并根据分辨率缩放
    u = NORMALIZE(u) * w;
    v = (1 - NORMALIZE(v)) * h;

    // 使用双线性插值计算纹理颜色
    int iu = (int)u, iv = (int)v;
    Vector3f ret_color = Vector3f::ZERO;
    float s = u - iu, t = v - iv;
    Vector3f color1 = (1 - s) * getColor(iu, iv + 1) + s * getColor(iu + 1, iv + 1);
    Vector3f color2 = (1 - s) * getColor(iu, iv) + s * getColor(iu + 1, iv);
    ret_color += (1 - t) * color2 + t * color1;
    return ret_color;
}

// 获取整数坐标点的颜色
Vector3f Texture::getColor(int x, int y) const {
    x = NORMALIZE2(x, w);
    y = NORMALIZE2(y, h);
    int index = (y * w + x) * channel;
    return Vector3f(pic[index + 0], pic[index + 1], pic[index + 2]) / 255.0;
}