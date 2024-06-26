#include "mesh.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>

void Mesh::computeNormal() {
    n.resize(t.size());
    for (int triId = 0; triId < (int)t.size(); ++triId) {
        TriangleIndex &triIndex = t[triId];
        Vector3f a = v[triIndex[1]] - v[triIndex[0]];
        Vector3f b = v[triIndex[2]] - v[triIndex[0]];
        b = Vector3f::cross(a, b);
        n[triId] = b / b.length();
    }
}
