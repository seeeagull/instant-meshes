#pragma once

#include <vector>

#ifndef INSTANTMESH_EXPORT
#  if defined(_MSC_VER) || defined(__CYGWIN__)
#    ifdef InstantMesh_EXPORT
#      define INSTANTMESH_EXPORT __declspec(dllexport)
#    else
#      define INSTANTMESH_EXPORT __declspec(dllimport)
#    endif
#  elif defined(__clang__) || defined(__GNUC__)
#    define INSTANTMESH_EXPORT __attribute__((visibility("default")))
#  endif
#endif


INSTANTMESH_EXPORT int runInstantMeshes(std::vector<std::vector<int>> &faces,
                            std::vector<std::vector<float>> &verts,
                            const std::vector<std::vector<int>> &features,
                            int argc, char **argv);
