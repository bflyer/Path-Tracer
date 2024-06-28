#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object3d.hpp"
#include "utils.hpp"
#include "plane.hpp"

#include <vecmath.h>
#include <cmath>
#include <iostream>
using namespace std;

class Triangle: public Object3D {

public:
	Triangle() = delete;

    // a b c are three vertex positions of the triangle
	Triangle( const Vector3f& a, const Vector3f& b, const Vector3f& c, Material* m)
		: Object3D(m), a(a), b(b), c(c), 
		an(Vector3f::ZERO), 
		bn(Vector3f::ZERO), 
		cn(Vector3f::ZERO),
		at(Vector2f::ZERO),
		bt(Vector2f::ZERO),
		ct(Vector2f::ZERO) {
		normal = Vector3f::cross((b - a), (c - a)).normalized();
		d = Vector3f::dot(normal, a);
        bound[0] = minE(minE(a, b), c);
        bound[1] = maxE(maxE(a, b), c);
        cen = (a + b + c) / 3;
        nSet = false;
        tSet = false;

		area = 0.5 * Vector3f::cross(b - a, c - a).length();
	}

	// // 法一：中心坐标法直接求交
	// bool intersect( const Ray& ray,  Hit& hit , float tmin) override {
	// 	Vector3f dir = ray.getDirection();
	// 	Vector3f e1 = a - b;
	// 	Vector3f e2 = a - c;
	// 	Vector3f s = a - ray.getOrigin();
	// 	Matrix3f m1(s, e1, e2, 1);
	// 	Matrix3f m2(dir, s, e2, 1);
	// 	Matrix3f m3(dir, e1, s, 1);
	// 	Matrix3f m4(dir, e1, e2, 1);
	// 	float d1 = m1.determinant();
	// 	float d2 = m2.determinant();
	// 	float d3 = m3.determinant();
	// 	float d4 = m4.determinant();
	// 	float t = d1/d4;
	// 	float b = d2/d4;
	// 	float c = d3/d4;
	// 	if (t <= hit.getT() && t >= tmin && b >= 0 && c >= 0 && b + c <= 1) {
	// 		hit.set(t, material, normal);
	// 		return true;
	// 	} else {
	// 		return false;
	// 	}
	// }

	// // 法二：先直接判断光线与三角形平面是否有交，然后判断交点是否在内部
	// bool intersect( const Ray& ray,  Hit& hit , float tmin) override {
	// 	float cos = Vector3f::dot(normal, ray.getDirection());
	// 	// 1. 判断光线与三角形平面是否有交
	// 	if (fabs(cos) < 1e-6) return false;
	// 	// 2. 找到交点
	// 	float d = Vector3f::dot(normal, a);
	// 	float t = (d - Vector3f::dot(normal, ray.getOrigin())) / cos;
	// 	if (t > hit.getT() || t < tmin) return false;
	// 	Vector3f point(ray.getOrigin() + ray.getDirection() * t);
	// 	// 3. 判断交点是否在三角形内部
	// 	if (!inTriangle(point)) return false;
	// 	else {
	// 		hit.set(t, material, normal);
	// 		return true;
	// 	}
	// }

	// 法三：叉乘判断法
	// Ref: ver.2020
	bool intersect(const Ray& ray, Hit& hit, float tmin) override {
        Vector3f o(ray.getOrigin()), dir(ray.getDirection());
		// 根据叉乘判断是否在三角形内有交
        Vector3f v0v1 = b - a;
        Vector3f v0v2 = c - a;
        Vector3f pvec = Vector3f::cross(dir, v0v2);
        float det = Vector3f::dot(v0v1, pvec);
        if (fabs(det) < tmin) return false;

		// 坐标变换
        float invDet = 1 / det;
        Vector3f tvec = o - a;
        float u = Vector3f::dot(tvec, pvec) * invDet;
        if (u < 0 || u > 1) return false;
        Vector3f qvec = Vector3f::cross(tvec, v0v1);
        float v = Vector3f::dot(dir, qvec) * invDet;
        if (v < 0 || u + v > 1) return false;
        double t = Vector3f::dot(v0v2, qvec) * invDet;
        if (t <= tmin || t > hit.getT()) return false;

		// 真正有交
        Vector3f p(o + dir * t);
        getUV(p, u, v);
        hit.set(t, material, getNorm(p), material->getColor(u, v), p);
        return true;
    }

	// Ref: ver.2020
	void setVNorm(const Vector3f& anorm, const Vector3f& bnorm,
                  const Vector3f& cnorm) {
        an = anorm;
        bn = bnorm;
        cn = cnorm;
        nSet = true;
    }

    void setVT(const Vector2f& _at, const Vector2f& _bt, const Vector2f& _ct) {
		at = _at;
        bt = _bt;
        ct = _ct;
        tSet = true;
    }

	double getArea() const override{
		return area;
	}

	Vector3f sample() const override {
		double r1 = RAND2;
		double r2 = (1 - r1) * RAND2;
		double r3 = 1 - r1 - r2;
		return r1 * a + r2 * b + r3 * c;
	}

	Vector3f min() const override { return bound[0]; }
    Vector3f max() const override { return bound[1]; }
    Vector3f center() const override { return cen; }
    vector<Object3D*> getFaces() override { return {(Object3D*)this}; }

	Vector3f normal;
	Vector3f a, b, c, cen;
    Vector2f at, bt, ct;
    Vector3f an, bn, cn;
    Vector3f bound[2];
	float d;
    bool nSet = false;
    bool tSet = false;

protected:
	// 判断一个点是否在三角形内部 
	bool inTriangle(const Vector3f& point) {
		return Vector3f::dot(Vector3f::cross(b - point, c - point), normal) >= -1e-6 &&
			   Vector3f::dot(Vector3f::cross(c - point, a - point), normal) >= -1e-6 &&
			   Vector3f::dot(Vector3f::cross(a - point, b - point), normal) >= -1e-6;
	}

	// Ref: ver.2020
	// 获取三角面片的法向量
	Vector3f getNorm(const Vector3f& p) {
        if (!nSet) return normal;
        Vector3f va = (a - p), vb = (b - p), vc = (c - p);
        float ra = Vector3f::cross(vb, vc).length(),
              rb = Vector3f::cross(vc, va).length(),
              rc = Vector3f::cross(va, vb).length();
        return (ra * an + rb * bn + rc * cn).normalized();
    }

	// Ref: ver.2020
	// 获取 u-v 坐标
    void getUV(const Vector3f& p, float& u, float& v) {
        if (!tSet) return;
        Vector3f va = (a - p), vb = (b - p), vc = (c - p);
        float ra = Vector3f::cross(vb, vc).length(),
              rb = Vector3f::cross(vc, va).length(),
              rc = Vector3f::cross(va, vb).length();
        Vector2f uv = (ra * at + rb * bt + rc * ct) / (ra + rb + rc);
        u = uv.x();
        v = uv.y();
    }

private:
	double area;
};

#endif 