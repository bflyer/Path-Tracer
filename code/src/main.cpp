#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <string>

#include "render.hpp"
#include "scene_parser.hpp"

using namespace std;

// int main(int argc, char *argv[]) {
//     for (int argNum = 1; argNum < argc; ++argNum) {
//         std::cout << "Argument " << argNum << " is: " << argv[argNum] << std::endl;
//     }

//     if (argc != 3) {
//         cout << "Usage: ./bin/PA1 <input scene file> <output bmp file>" << endl;
//         return 1;
//     }
//     string inputFile = argv[1];
//     string outputFile = argv[2];  // only bmp is allowed.
//     // Step 1: Parse the scene using SceneParser.
//     SceneParser sceneParser = SceneParser(inputFile.c_str());        // 获取场景解析器
//     Camera* camera = sceneParser.getCamera();                        // 获取相机
//     Group* baseGroup = sceneParser.getGroup();                       // 获取对象列表
//     Image outImg = Image(camera->getWidth(), camera->getHeight());    // 准备图像以供渲染

//     // Step 2: Loop over each pixel in the image, shooting a ray
//     // through that pixel and finding its intersection with
//     // the scene.  
//     for (int x = 0; x < camera->getWidth(); x++) {
//         for (int y = 0; y < camera->getHeight(); y++) {
//             // A. 计算当前像素 (x, y) 处相机出射光线 camRay
//             Ray camRay = camera->generateRay(Vector2f(x, y));
//             // 注意这里需要对每个像素用一个新的 hit ，避免之前的影响
//             Hit hit;
//             // B. 判断 camRay 是否和场景有交点，并返回最近交点的数据，存储在 hit 中
//             bool isIntersect = baseGroup->intersect(camRay, hit, 0);
//             if (isIntersect) {
//                 Vector3f finalColor = Vector3f::ZERO;
//                 // C.2 找到交点之后，累加来自所有光源的光强影响
//                 // Step 3: Write the color at the intersection to that
//                 // pixel in your output image.
//                 for (int li = 0; li < sceneParser.getNumLights(); li++) {
//                     Light* light = sceneParser.getLight(li);
//                     Vector3f L;
//                     Vector3f lightColor;
//                     // a. 获得光源
//                     light->getIllumination(camRay.pointAtParameter(hit.getT()), L, lightColor);
//                     // b. 计算局部光强
//                     finalColor += hit.getMaterial()->Shade(camRay, hit, L, lightColor);
//                 }
//                 outImg.SetPixel(x, y, finalColor);
//             } else {
//                 // C.1 不存在交点，返回背景色
//                 outImg.SetPixel(x, y, sceneParser.getBackgroundColor());
//             }
//         }
//     }
//     cout << "Hello! Computer Graphics!" << endl;
    
//     // Step 4: Save the output image to a file.
//     outImg.SaveBMP(outputFile.c_str());
//     return 0;
// }

// Ref: ver.2020
int main(int argc, char *argv[]) {
    for (int argNum = 1; argNum < argc; ++argNum) {
        std::cout << "Argument " << argNum << " is: " << argv[argNum] << std::endl;
    }
    // if (argc < 4) {
    //     std::cout << "Usage: ./bin/PA1 <input scene file> <output bmp file> !!!"
    //                 "<method> <spp>"
    //               << endl;
    //     return 1;
    // }

    SceneParser sceneParser(argv[1]);

    if (!strcmp(argv[3], "rc") || !strcmp(argv[3], "pt") || !strcmp(argv[3], "wt")) {
        int samps = atoi(argv[4]);
        PathTracer pt(sceneParser, samps, argv[2], argv[3]);
        pt.render();
    } else {
        cout << "Unknown method: " << argv[3] << endl;
        return 1;
    }
    return 0;
}