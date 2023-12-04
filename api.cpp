#include "src/batch.h"
#include "src/viewer.h"
#include "src/serializer.h"
#include "api.h"
#include <thread>
#include <cstdlib>

/* Force usage of discrete GPU on laptops */
NANOGUI_FORCE_DISCRETE_GPU();

int runInstantMeshes(const std::vector<std::vector<int>> &faces,
                     const std::vector<std::vector<float>> &verts,
                     const std::vector<std::vector<int>> &features,
                     int argc, char **argv) {
    int nprocs = -1;
    std::vector<std::string> args;
    bool extrinsic = true, dominant = false, align_to_boundaries = false;
    bool fullscreen = false, help = false, deterministic = false, compat = false;
    int rosy = 4, posy = 4, face_count = -1, vertex_count = -1;
    uint32_t knn_points = 10, smooth_iter = 2;
    Float crease_angle = -1, scale = -1;
    std::string batchOutput;
    #if defined(__APPLE__)
        bool launched_from_finder = false;
    #endif

    try {
        for (int i=1; i<argc; ++i) {
            if (strcmp("--deterministic", argv[i]) == 0 || strcmp("-d", argv[i]) == 0) {
                deterministic = true;
            } else if (strcmp("--intrinsic", argv[i]) == 0 || strcmp("-i", argv[i]) == 0) {
                extrinsic = false;
            } else if (strcmp("--boundaries", argv[i]) == 0 || strcmp("-b", argv[i]) == 0) {
                align_to_boundaries = true;
            } else if (strcmp("--smooth", argv[i]) == 0 || strcmp("-S", argv[i]) == 0) {
                smooth_iter = str_to_uint32_t(argv[++i]);
            } else if (strcmp("--knn", argv[i]) == 0 || strcmp("-k", argv[i]) == 0) {
                knn_points = str_to_uint32_t(argv[++i]);
            } else if (strcmp("--crease", argv[i]) == 0 || strcmp("-c", argv[i]) == 0) {
                crease_angle = str_to_float(argv[++i]);
            } else if (strcmp("--rosy", argv[i]) == 0 || strcmp("-r", argv[i]) == 0) {
                rosy = str_to_int32_t(argv[++i]);
            } else if (strcmp("--posy", argv[i]) == 0 || strcmp("-p", argv[i]) == 0) {
                posy = str_to_int32_t(argv[++i]);
                if (posy == 6)
                    posy = 3;
            } else if (strcmp("--scale", argv[i]) == 0 || strcmp("-s", argv[i]) == 0) {
                scale = str_to_float(argv[++i]);
            } else if (strcmp("--faces", argv[i]) == 0 || strcmp("-f", argv[i]) == 0) {
                face_count = str_to_int32_t(argv[++i]);
            } else if (strcmp("--vertices", argv[i]) == 0 || strcmp("-v", argv[i]) == 0) {
                vertex_count = str_to_int32_t(argv[++i]);
            } else if (strcmp("--output", argv[i]) == 0 || strcmp("-o", argv[i]) == 0) {
                batchOutput = argv[++i];
            } else if (strcmp("--dominant", argv[i]) == 0 || strcmp("-D", argv[i]) == 0) {
                dominant = true;
            } else if (strcmp("--compat", argv[i]) == 0 || strcmp("-C", argv[i]) == 0) {
                compat = true;
#if defined(__APPLE__)
            } else if (strncmp("-psn", argv[i], 4) == 0) {
                launched_from_finder = true;
#endif
            } else {
                args.push_back(argv[i]);
            }
        }
    } catch (const std::exception &e) {
        cout << "Error: " << e.what() << endl;
        help = true;
    }

    if ((posy != 3 && posy != 4) || (rosy != 2 && rosy != 4 && rosy != 6)) {
        cerr << "Error: Invalid symmetry type!" << endl;
        help  = true;
    }

    tbb::task_scheduler_init init(nprocs == -1 ? tbb::task_scheduler_init::automatic : nprocs);

    if (!batchOutput.empty() && args.size() == 1) {
        try {
            batch_process(args[0], batchOutput, rosy, posy, scale, face_count,
                          vertex_count, crease_angle, extrinsic,
                          align_to_boundaries, smooth_iter, knn_points,
                          !dominant, deterministic, faces, verts, features);
            return 0;
        } catch (const std::exception &e) {
            cerr << "Caught runtime error : " << e.what() << endl;
            return -1;
        }
    }

    return EXIT_SUCCESS;
}