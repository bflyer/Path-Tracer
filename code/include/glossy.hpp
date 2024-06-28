#include <vecmath.h>
#include <math.h>

#define minf(a, b) (((a) < (b)) ? (a) : (b))
#define maxf(a, b) (((a) > (b)) ? (a) : (b))
#define EPSILON 1e-5f  // 用于分母中以确保不会发生除以零的操作

// GGX分布函数 (D)
float D_GGX(float NoH, float alpha) {
    float alphaSqr = alpha * alpha;
    float NoH2 = NoH * NoH;
    float denom = NoH2 * alphaSqr + (1.0 - NoH2);
    return alphaSqr / (M_PI * denom * denom);
}

// 几何项 (G)
float G_Smith(float NoV, float NoL, float alpha) {
    float a = alpha;
    float GGXV = NoL * (NoV * (1.0 - a) + a);
    float GGXL = NoV * (NoL * (1.0 - a) + a);
    return 0.5 / (GGXV + GGXL);
}

// Fresnel项 (F)，这里使用Schlick近似
Vector3f F_Schlick(float VoH, Vector3f F0) {
    float f = pow(1.0 - VoH, 5.0);
    return F0 + (Vector3f(1.0) - F0) * f;
}

Vector3f mix(Vector3f a, Vector3f b, float t) {
    return t * a + (1 - t) * b;
}

/*
    最终的 BRDF 计算
    N: 法线
    V: 视线方向
    L: 光线方向
    albedo: 表面的基色或漫反射颜色
    metallic: 金属度参数，取值范围通常是[0, 1]
            (0 表示完全的非金属材料（如塑料、木头），
            1 表示完全的金属材料（如金、银）)
    roughness: 粗糙度参数，同样取值范围是[0, 1]，
            描述了表面微观结构的平滑程度。
            0 表示非常光滑的表面，反射高光会非常集中和锐利；
            而 1 表示非常粗糙的表面，高光会分散并变得更加模糊。
    常见例子：
        albedo:
            对于非金属材质，如草地，可能有(0.25, 0.5, 0.1)（暗绿色）的albedo。
            金属材质如铜，可能会用(0.75, 0.45, 0.2)（带有红色和橙色调）
            的 albedo，并且其金属度会设置得较高。
        metallic:
            非金属材质如布料或皮肤，metallic 值通常设为 0。
            金属材质如铁或银，metallic 值则设为 1。
        roughness:
            镜子表面的 roughness 接近0，因为它几乎完全光滑。
            橡胶或未打磨的金属表面可能有较高的 roughness 值，
            比如 0.7 到 0.9，这会让它们的高光更加散开和模糊。
*/
Vector3f BRDF_GGX(Vector3f N, Vector3f V, Vector3f L, Vector3f albedo, float metallic, float roughness) {
    // 视线方向和光线方向到法线的夹角余弦值
    float NoV = maxf(Vector3f::dot(N, V), 0.0);
    float NoL = maxf(Vector3f::dot(N, L), 0.0);
    
    // 半向量H
    Vector3f H = (V + L).normalized();
    float NoH = maxf(Vector3f::dot(N, H), 0.0);
    float VoH = maxf(Vector3f::dot(V, H), 0.0);
    
    // F0，金属的基色反射率
    Vector3f F0 = Vector3f(0.04); // 默认为非金属的菲涅尔基础值（约等于塑料的）
    F0 = mix(F0, albedo, metallic); // 根据metallic混合
    
    // 计算各项
    float D = D_GGX(NoH, roughness);
    float G = G_Smith(NoV, NoL, roughness);
    Vector3f F = F_Schlick(VoH, F0);
    
    // 微面元的双向反射分布函数
    Vector3f specular = (D * G * F) / (4.0 * NoV * NoL + EPSILON); // EPSILON是一个小常数，用于防止除以零
    
    // Diffuse部分，对于金属，diffuse贡献为0
    Vector3f diffuse = albedo * (1.0 - metallic) * (1.0 / M_PI);
    
    return diffuse + specular;
}