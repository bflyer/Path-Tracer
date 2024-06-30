#ifndef CURVE_HPP
#define CURVE_HPP

#include "object3d.hpp"
#include <utility>
#include <vecmath.h>
#include <vector>

#include <algorithm>

// Ref: ver.2023
// TODO (PA2): Implement Bernstein class to compute spline basis function.
//       You may refer to the python-script for implementation.

// The CurvePoint object stores information about a point on a curve
// after it has been tesselated: the vertex (V) and the tangent (T)
// It is the responsiblility of functions that create these objects to fill in all the data.
struct CurvePoint {
    Vector3f V; // Vertex
    Vector3f T; // Tangent  (unit)
};

class Curve : public Object3D {
public:
    explicit Curve(std::vector<Vector3f> points)
        : controls(std::move(points)) {
        ymin = INF;
        ymax = -INF;
        radius = 0;
        for (auto pt : controls) {
            ymin = std::min(pt.y(), ymin);
            ymax = std::max(pt.y(), ymax);
            radius = std::max(radius, fabs(pt.x()));
            radius = std::max(radius, fabs(pt.z()));
        }
    }

    bool intersect(const Ray& r, Hit& h, float tmin) override { return false; }

    std::vector<Vector3f>& getControls() { return controls; }

    std::pair<double,double> get_valid_range(){ return std::make_pair(t[k],t[n]); }
    void evaluate(double mu){
        // 使用 upper_bound 找到第一个大于 mu 的元素在 t 中的位置，
        // 然后减 1 得到不大于 mu 的最大元素下标
        int bpos = upper_bound(t.begin(), t.end(), mu) - t.begin() - 1;
        s = std::vector<double> (k + 2, 0), ds = std::vector<double> (k + 1, 1);
        s[k] = 1;  // 递推基
        for (int p = 1; p < k + 1; p++) {
            // 内层循环遍历从 k - p 到 k，结合 bpos 使用以保持索引合法
            for (int ii = k - p; ii < k + 1; ii++) {
                // 调整索引 i 以适应原始数据 tpad
                int i = ii + bpos - k;
                // 计算 w1, dw1, w2, dw2，这些是基于 mu、tpad 的差值及比值的权重和导数
                double w1, dw1, w2, dw2;
                if (tpad[i + p] == tpad[i]) {
                    w1 = mu;
                    dw1 = 1;
                } else {
                    w1 = (mu - tpad[i]) / (tpad[i + p] - tpad[i]);
                    dw1 = 1.0 / (tpad[i + p] - tpad[i]);
                }
                if (tpad[i + p + 1] == tpad[i + 1]) {
                    w2 = 1 - mu;
                    dw2 = -1;
                } else {
                    w2 = (tpad[i + p + 1] - mu) / (tpad[i + p + 1] - tpad[i + 1]);
                    dw2 = -1 / (tpad[i + p + 1] - tpad[i + 1]);
                }
                // 当 p 等于 k 时，直接计算 ds 的值
                if (p == k)
                    ds[ii] = (dw1 * s[ii] + dw2 * s[ii + 1]) * p;
                // 基于 w1, w2 以及前一状态的 s 值更新 s 的值
                s[ii] = w1 * s[ii] + w2 * s[ii + 1];
            }
        }
        // 删除 s 的最后一个元素，因为其不参与最终计算
        s.pop_back();
        // 计算左右两侧需要保留的元素数量，用于后续裁剪操作
        lsk = bpos-k; 
        int rsk = n-bpos - 1;
        // 如果 lsk 小于 0，说明左侧需要移除一些元素以保持合法索引
        if (lsk < 0) {
            s.erase(s.begin(),s.begin()-lsk);
            ds.erase(ds.begin(),ds.begin()-lsk);
            lsk = 0;
        }
        // 同理，如果 rsk 小于 0，右侧也需要进行相应元素的移除
        if (rsk < 0) {
            s.erase(s.end()+rsk,s.end());
            ds.erase(ds.end()+rsk,ds.end());
        }
    }

    void pad() {
        int tSize = t.size();
        tpad.resize(tSize + k);
        for (int i = 0; i < tSize + k; i++) 
            tpad[i] = i < tSize ? t[i] : t.back();
    }
    virtual void discretize(int resolution, std::vector<CurvePoint>& data) = 0;

    std::vector<Vector3f> controls;
    int n, k; // n: number of control points, k: number of degree
    std::vector<double> t, tpad; // t: knot vector
    std::vector<double> s, ds;
    float ymin, ymax, radius;
    int lsk;
    double range[2];
};

class BezierCurve : public Curve {
public:
    explicit BezierCurve(const std::vector<Vector3f>& points)
        : Curve(points)
    {
        if (points.size() < 4 || points.size() % 3 != 1) {
            printf("Number of control points of BezierCurve must be 3n+1!\n");
            exit(0);
        }
        n = controls.size();
        k = n - 1;
        range[0] = 0;
        range[1] = 1;
        t.resize(2 * n);
        for (int i = 0; i < n; i++) {
            t[i] = 0;
            t[i + n] = 1;
        }
        pad();
    }

    // 初始化组合数（填充边界）
    vector<vector<int>> initCombinationTable(int n) {
        vector<vector<int>> Comb(n + 1, vector<int>(n + 1, 0));
        for (int i = 0; i <= n; ++i) {
            Comb[i][0] = 1; // C(i, 0) 总是1
            Comb[i][i] = 1; // C(i, i) 总是1
        }
        return Comb;
    }

    // 动态规划计算组合数
    void computeCombinations() {
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j < i; ++j) {
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
            }
        }
    }

    // 生成顶点系数
    void computeBpre(float t) {
        for (int i = 0; i <= n-1; i++) {
            Bpre.push_back(C[n-1][i] * pow(t, i) * pow(1-t, n-1-i));
        }
    }

    // 生成顶点向量（Pi * Bi 加和）
    Vector3f fV(float t) {
        auto sum = Vector3f::ZERO;
        // 首尾直接计算
        sum += controls[0] * C[n][0] * pow(1-t, n);
        // 中间利用 B(n-1, i-1)(t) 和 B(n-1, i)(t) 计算
        for (int i = 1; i < n; ++i) {
            auto Pi = controls[i];
            // B = (1-t) * B(n-1, i)(t) + t * B(n-1, i-1)(t)
            auto Bint = (1-t) * Bpre[i] + t * Bpre[i-1];
            sum += Pi * Bint;
        }
        sum += controls[n] * C[n][n] * pow(t, n);
        return sum;
    }

    // 生成切向量
    Vector3f fT(float t) {
        auto sum = Vector3f::ZERO;
        // 首尾直接计算
        sum += controls[0] * (-float(n) * pow(1 - t, n - 1));
        // 中间利用 B(n-1, i-1)(t) 和 B(n-1, i)(t) 计算
        for (int i = 1; i < n; ++i) {
            auto Pi = controls[i];
            // B'(n, i)(t) = n * [B(n-1, i-1)(t) - B(n-1, i)(t)]
            auto B_prime_int = n * (Bpre[i-1] - Bpre[i]);
            sum += Pi * B_prime_int;
        }
        sum += controls[n] * float(n) * pow(t, n - 1);
        return sum.normalized();
    }

    void discretize(int resolution, std::vector<CurvePoint>& data) override {
        data.clear();
        // TODO (PA2): fill in data vector
        // 对于三次 Bezier 分段连续曲线，resolution 指每一段连续区间 [t_{3i},t_{3i+3}] 离散化出的子区间个数
        // 也就是在 t = 3i/n 和 t = 3i/n + 3/n 之间，要取样 resolution 次 
        // 第 j 个取样点，对应的 t 就应该是 3i/n + 3j/(resolution * n)
        auto delta = 1.0f / float(resolution);
        // 注意 n 要从 0 到 resolution 而不是 resolution - 1
        for (int i = 0; i <= resolution; ++i) {
            Bpre.clear();
            computeBpre(i * delta);
            data.push_back(CurvePoint{fV(i * delta), fT(i * delta)});
        }
    }

    std::vector<vector<int>> C;  // 组合数
    std::vector<float> Bpre;       // B(i, n-1)(t)
};

class BsplineCurve : public Curve {
public:
    BsplineCurve(const std::vector<Vector3f>& points)
        : Curve(points) {
        if (points.size() < 4) {
            printf("Number of control points of BspineCurve must be more than 4!\n");
            exit(0);
        }
        n = controls.size();
        k = 3;
        t.resize(n + k + 1);
        for (int i = 0; i < n + k + 1; ++i)
            t[i] = (double)i / (n + k);
        pad();
        range[0] = t[k];
        range[1] = t[n];
    }

    // 得到均匀 B 样条的节点向量
    float ti(int i) {
        assert(i >= 0 and i <= n + k + 1);
        return float(i) / float(k + n + 1);
    }

    // 自底向上动态规划实现的B样条函数
    std::vector<std::vector<float>> calculateBSplineBasisFunctions(int n, int k, float t) {
        std::vector<std::vector<float>> dp(n + k + 2, std::vector<float>(k + 1)); // 初始化dp数组

        // 处理基本情况 p = 0
        for (int i = 0; i <= n + k + 1; ++i) {
            float start = ti(i), end = ti(i + 1);
            if (t >= start && t < end) {
                dp[i][0] = 1.0f;
            } else {
                dp[i][0] = 0.0f;
            }
        }

        // 自底向上计算更高阶的B样条基函数
        for (int p = 1; p <= k; ++p) {
            for (int i = 0; i <= n + k - p; ++i) {
                float denom1 = ti(i + p) - ti(i);
                float denom2 = ti(i + p + 1) - ti(i + 1);
                if (denom1 == 0 || denom2 == 0) continue; // 防止除以零错误
                
                dp[i][p] = 
                    ((t - ti(i)) / denom1) * dp[i][p - 1] +
                    ((ti(i + p + 1) - t) / denom2) * dp[i + 1][p - 1];
            }
        }

        // 返回最高阶的B样条基函数（通常只返回所需的那一阶）
        return dp;
    }

    // 计算基函数
    float B(int i, int p, float t) {
        return DP[i][p];
    }

    // 计算积函数的导数
    float Bd(int i, int p, float t) {
        return B(i, p - 1, t) - B(i + 1, p - 1, t);
    }

    // 生成顶点向量
    Vector3f fV(float t) {
        auto sum = Vector3f::ZERO;
        for (int i = 0; i <= n; ++i) {
            auto Pi = controls[i];
            auto Bikt = B(i, k, t);
            sum += Pi * Bikt;
        }
        return sum;
    }

    // 生成切向量
    Vector3f fT(float t) {
        auto sum = Vector3f::ZERO;
        for (int i = 0; i <= n; ++i) {
            auto Pi = controls[i];
            auto Bikt = Bd(i, k, t);
            sum += Pi * Bikt;
        }

        return sum.normalized();
    }

    void discretize(int resolution, std::vector<CurvePoint> &data) override {
        data.clear();
        // TODO (PA2): fill in data vector
        n = controls.size() - 1;
        k = 3;
        for (int i = k; i <= n; ++i) {
            auto delta = (ti(i + 1) - ti(i)) / resolution;
            float lo = ti(i);

            for (int j = 0; j <= resolution; ++j) {
                DP = calculateBSplineBasisFunctions(n, k, j * delta + lo);
                data.push_back(CurvePoint{fV(j * delta + lo), fT(j * delta + lo)});
            }
        }
        // // 补上尾巴
        // float hi = ti(n + 1);
        // data.push_back(CurvePoint{fV(hi), fT(hi)});
    }
protected:
    int n;  // 控制点数量
    int k;  // 阶数
    std::vector<std::vector<float>> DP;
};

#endif // CURVE_HPP