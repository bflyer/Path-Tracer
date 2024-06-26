#ifndef CAMERA_H
#define CAMERA_H

#include "ray.hpp"
#include <vecmath.h>
#include <float.h>
#include <cmath>


class Camera {
public:
    Camera(const Vector3f &center, const Vector3f &direction, const Vector3f &up, int imgW, int imgH) {
        this->center = center;                                                      // 光心
        this->direction = direction.normalized();                                   // 相机视线方向
        this->horizontal = Vector3f::cross(this->direction, up).normalized();                    // 水平轴
        this->up = Vector3f::cross(this->horizontal, this->direction).normalized();              // up 方向
        this->width = imgW;                                                         // 成像屏幕宽度
        this->height = imgH;                                                        // 成像屏幕高度
    }

    // Generate rays for each screen-space coordinate
    virtual Ray generateRay(const Vector2f &point) = 0;
    virtual ~Camera() = default;

    int getWidth() const { return width; }
    int getHeight() const { return height; }

    void setCenter(const Vector3f& pos) {
        this->center = pos;
    }
    Vector3f getCenter() const {
        return this->center;
    }

    void setRotation(const Matrix3f& mat) {
        this->horizontal = mat.getCol(0);
        this->up = -mat.getCol(1);
        this->direction = mat.getCol(2);
    }
    Matrix3f getRotation() const {
        return Matrix3f(this->horizontal, -this->up, this->direction);
    }

    virtual void resize(int w, int h) {
        width = w; height = h;
    }

protected:
    // Extrinsic parameters
    Vector3f center;
    Vector3f direction;
    Vector3f up;
    Vector3f horizontal;
    // Intrinsic parameters
    int width;
    int height;
};

// You can add new functions or variables whenever needed.
class PerspectiveCamera : public Camera {
public:
    PerspectiveCamera(const Vector3f &center, const Vector3f &direction,
            const Vector3f &up, int imgW, int imgH, float angle,
            float f = 20.0f, float aperture = 1.0f)
            : Camera(center, direction, up, imgW, imgH),
            focalLength(f),
            aperture(aperture) {
        // angle is fovy in radian.
        fovyd = angle / 3.1415 * 180.0;
        fx = fy = (float) height / (2 * tanf(angle / 2));
        cx = width / 2.0f;                                   // 光心坐标
        cy = height / 2.0f;
    }

    void resize(int w, int h) override {
        fx *= (float) h / height;
        fy = fx;
        Camera::resize(w, h);
        cx = width / 2.0f;
        cy = height / 2.0f;
    }

    // TODO：这里换了相机，要小心
    // Ray generateRay(const Vector2f &point) override {
    //     Vector3f d_rc = Vector3f((point[0] - cx) / fx, (cy - point[1]) / fy, 1);
    //     c2w = Matrix3f(horizontal, -up, direction, 1);  // 相机坐标系到世界坐标系的变换矩阵
    //     return Ray(center, c2w * d_rc);
    // }
    Ray generateRay(const Vector2f &point) override {
        float csx = (point.x() - cx) / fx;
        float csy = (point.y() - cy) / fy;
        Vector3f dir(csx, -csy, 1.0f);
        Matrix3f R(this->horizontal, -this->up, this->direction);
        dir = R * dir;
        dir = dir / dir.length();
        Ray ray(this->center, dir);
        return ray;
    }

private:
    // Perspective intrinsics
    float fx;
    float fy;
    float cx;
    float cy;
    float fovyd;
    float aperture;         // 光圈
    float focalLength;      // 焦距
    Matrix3f c2w;
};
#endif //CAMERA_H