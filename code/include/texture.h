#ifndef TEXTURE_H
#define TEXTURE_H
// Ref: ver.2023
#include <string>

#include <vecmath.h>
using std::string;
struct Texture {  // 纹理
    unsigned char *pic;
    int w, h, channel;           // 宽、高、通道数
    Texture(const char *textureFile);
    
    // 由 (u, v) 坐标获取颜色
    Vector3f getColor(float u, float v) const;
    Vector3f getColor(int u, int v) const;

    // 获取扰动
    float getDisturb(float u, float v, Vector2f &grad) const;
    inline int getIdx(float u, float v) const;
    inline float getGray(int idx) const;
};

#endif  // !TEXTURE_H