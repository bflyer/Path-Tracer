// // Path Tracing（递推版）
// // Ref: smallpt
// static Vector3f ptColor(Ray ray, const SceneParser& sceneParser, int depth = 0, int E = 1) {
//     Group* group = sceneParser.getGroup();
//     Vector3f color = Vector3f::ZERO;
//     Vector3f opacity = Vector3f(1, 1, 1);
//     Hit hit;
//     Ray r = ray;

//     while (true) {
//         // 1. 求交：如果没有交点，直接返回背景色
//         // TODO: 返回背景色还是返回黑色呢？
//         if (!group->intersect(r, hit, TMIN)) {
//             return color;
//         }
        
//         // 交点坐标
//         Vector3f hitPos(r.getOrigin() + r.getDirection() * hit.getT());
//         Material* material = hit.getMaterial();          // the hit object
//         Vector3f f(hit.getColor());         // BRDF
//         Vector3f n(hit.getNormal());
//         Vector3f nl = Vector3f::dot(n, r.getDirection()) < 0 ? n : -n;  // ensure the normal points outward
        
//         // TODO：丢弃时选择 0 还是选择材质的 emission
//         // 2. R.R.(Russian Roulette)
//         if (++depth > RR_DEPTH) {
//             if (RAND2 < RR_PROBABILITY)
//                 f = f * (1.0 / RR_PROBABILITY);
//             else 
//                 return color;
//         }

//         color = color + opacity * material->getEmission();
//         opacity = opacity * f;

//         float type = RAND2;
//         // 3. 单独处理理想漫反射和理想镜面反射
//         // Ideal DIFFUSE reflection(理想漫反射)
//         if (material->getType().x() == 1){             
//             // 记极角(polar angle) 为 theta，记方位角(Azimuthal angle) 为 phi    
//             double phi = 2 * M_PI * RAND2;  
//             double sinTheta2 = RAND2;              // sin(theta) ^ 2
//             double sinTheta = sqrt(sinTheta2);     // sin(theta)
            
//             // 构建局部坐标系
//             Vector3f w = nl;  // 表面的单位法线
//             // 通过判断法线向量 w 的 x 分量绝对值是否大于 0.1 来决定
//             // 构造正交基的方式，以避免因法线几乎平行于 x 轴而导致的
//             // 除以零问题。
//             Vector3f u = (Vector3f::cross((fabs(w.x()) > .1 ? Vector3f(0, 1, 0) : Vector3f(1, 0, 0)), w)).normalized();
//             Vector3f v = Vector3f::cross(w, u);
//             // 利用随机数 phi 和 sinTheta，以及正交基 u 和 v，
//             // 通过球坐标系转换为笛卡尔坐标系的方式计算出新的反射方向向量 d
//             Vector3f d = (u * cos(phi) * sinTheta + v * sin(phi) * sinTheta + w * sqrt(1 - sinTheta2)).normalized();

//             r = Ray(hitPos, d);
//             continue;
//         }

//         // Ideal SPECULAR reflection(理想镜面反射)
//         else if (material->getType().y() == 1) {
//             r = Ray(hitPos, r.getDirection() - n * 2 * Vector3f::dot(r.getDirection(), n));
//             continue;
//         }
            
//         // Ideal dielectric REFRACTION(理想介质折射)
//         // 反射光线初始化，直接套用反射公式
//         Ray reflRay(hitPos, ray.getDirection() - n * 2 * Vector3f::dot(ray.getDirection(), n)); 
//         // into = true: 光线从外而内；into = false：光线从内而外    
//         bool into = Vector3f::dot(n, nl) > 0;                
//         // nc: 外部介质折射率（如空气）
//         // nt: 内部介质的折射率（如玻璃）
//         double nc = 1, nt = material->getRefractRate();
//         // nnt: 根据光线进出方向计算相对折射率
//         double nnt = into ? nc/nt : nt/nc;
//         // ddn: 入射角余弦值
//         double ddn = Vector3f::dot(ray.getDirection(), nl);
//         double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
//         // Total internal reflection(全内反射检查)
//         if (cos2t < 0) {
//             r = reflRay;
//             continue;
//         }
            
//         // 计算折射方向
//         Vector3f tdir = (ray.getDirection() * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalized();
//         // 计算菲涅尔反射和折射系数
//         // R0: 菲涅耳反射系数，最终反射率
//         // Re: 最终折射率；Tr: 最终反射率
//         // c: 辅助变量，用于确定光线能量在反射和折射间的分配
//         double a = nt - nc, b = nt + nc, R0 = a * a / (b*b), c = 1 - (into ? -ddn : Vector3f::dot(tdir, n));
//         double Re = R0 + (1 - R0) *c*c*c*c*c;
//         double Tr = 1 - Re;
//         // P: 用于决定追踪反射光线还是折射光线
//         // RP: 追踪反射光线的权重；TP: 追踪折射光线的权重
//         double P = .25 + .5*Re, RP = Re/P, TP = Tr / (1-P);
//         // R.R.
//         // 递归深度大于阈值时使用俄罗斯轮盘赌
//         if (RAND2 < P) {
//             opacity = opacity * RP;
//             r = reflRay;
//         } else {
//             opacity = opacity * TP;
//             r = Ray(hitPos, tdir);
//         }
//     }
// }