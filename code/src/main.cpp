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

// Ref: ver.2020
int main(int argc, char *argv[]) {
    for (int argNum = 1; argNum < argc; ++argNum) {
        std::cout << "Argument " << argNum << " is: " << argv[argNum] << std::endl;
    }
    if (argc < 4) {
        std::cout << "Usage: ./bin/PA1 <input scene file> <output bmp file> !!!"
                    "<method> <spp>"
                  << endl;
        return 1;
    }

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