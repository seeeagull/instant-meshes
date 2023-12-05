/*
    batch.cpp -- command line interface to Instant Meshes

    This file is part of the implementation of

        Instant Field-Aligned Meshes
        Wenzel Jakob, Daniele Panozzo, Marco Tarini, and Olga Sorkine-Hornung
        In ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2015)

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE.txt file.
*/

#include "batch.h"
#include "meshio.h"
#include "dedge.h"
#include "subdivide.h"
#include "meshstats.h"
#include "hierarchy.h"
#include "field.h"
#include "normal.h"
#include "extract.h"
#include "bvh.h"

void batch_process(const std::string &input, const std::string &output,
                   int rosy, int posy, Float scale, int face_count,
                   int vertex_count, Float creaseAngle, bool extrinsic,
                   bool align_to_boundaries, int smooth_iter, int knn_points,
                   bool pure_quad, bool deterministic,
                   std::vector<std::vector<int>> &faces,
                   std::vector<std::vector<float>> &verts,
                   const std::vector<std::vector<int>> &features) {
    cout << endl;
    cout << "Running in batch mode:" << endl;
    cout << "   Input file             = " << input << endl;
    cout << "   Output file            = " << output << endl;
    cout << "   Rotation symmetry type = " << rosy << endl;
    cout << "   Position symmetry type = " << (posy==3?6:posy) << endl;
    cout << "   Crease angle threshold = ";
    if (creaseAngle > 0)
        cout << creaseAngle << endl;
    else
        cout << "disabled" << endl;
    cout << "   Extrinsic mode         = " << (extrinsic ? "enabled" : "disabled") << endl;
    cout << "   Align to boundaries    = " << (align_to_boundaries ? "yes" : "no") << endl;
    cout << "   kNN points             = " << knn_points << " (only applies to point clouds)"<< endl;
    cout << "   Fully deterministic    = " << (deterministic ? "yes" : "no") << endl;
    if (posy == 4)
        cout << "   Output mode            = " << (pure_quad ? "pure quad mesh" : "quad-dominant mesh") << endl;
    cout << endl;

    MatrixXu F;
    MatrixXf V, N;
    VectorXf A;
    std::set<uint32_t> crease_in, crease_out;
    BVH *bvh = nullptr;
    AdjacencyMatrix adj = nullptr;

    if (input.size() > 0) {
        /* Load the input mesh */
        load_mesh_or_pointcloud(input, F, V, N);
    } else {
        int face_num = faces.size(), vert_num = verts.size();
        F.conservativeResize(3, face_num);
        V.conservativeResize(3, vert_num);
        for (int i = 0; i < face_num; ++i) {
            for (int j = 0; j < 3; ++j)
                F(j, i) = faces[i][j];
        }
        for (int i = 0; i < vert_num; ++i) {
            for (int j = 0; j < 3; ++j)
                V(j, i) = verts[i][j];
        }
    }

    bool pointcloud = F.size() == 0;

    Timer<> timer;
    MeshStats stats = compute_mesh_stats(F, V, deterministic);

    if (pointcloud) {
        bvh = new BVH(&F, &V, &N, stats.mAABB);
        bvh->build();
        adj = generate_adjacency_matrix_pointcloud(V, N, bvh, stats, knn_points, deterministic);
        A.resize(V.cols());
        A.setConstant(1.0f);
    }

    if (scale < 0 && vertex_count < 0 && face_count < 0) {
        cout << "No target vertex count/face count/scale argument provided. "
                "Setting to the default of 1/16 * input vertex count." << endl;
        vertex_count = V.cols() / 16;
    }

    if (scale > 0) {
        Float face_area = posy == 4 ? (scale*scale) : (std::sqrt(3.f)/4.f*scale*scale);
        face_count = stats.mSurfaceArea / face_area;
        vertex_count = posy == 4 ? face_count : (face_count / 2);
    } else if (face_count > 0) {
        Float face_area = stats.mSurfaceArea / face_count;
        vertex_count = posy == 4 ? face_count : (face_count / 2);
        scale = posy == 4 ? std::sqrt(face_area) : (2*std::sqrt(face_area * std::sqrt(1.f/3.f)));
    } else if (vertex_count > 0) {
        face_count = posy == 4 ? vertex_count : (vertex_count * 2);
        Float face_area = stats.mSurfaceArea / face_count;
        scale = posy == 4 ? std::sqrt(face_area) : (2*std::sqrt(face_area * std::sqrt(1.f/3.f)));
    }

    cout << "Output mesh goals (approximate)" << endl;
    cout << "   Vertex count           = " << vertex_count << endl;
    cout << "   Face count             = " << face_count << endl;
    cout << "   Edge length            = " << scale << endl;

    MultiResolutionHierarchy mRes;
    MatrixXu Feat(3, F.cols());
    for (int i = 0; i < features.size(); ++i)
        for (int j = 0; j < 3; ++j) {
            Feat(j, i) = features[i][j];
        }

    if (!pointcloud) {
        /* Subdivide the mesh if necessary */
        VectorXu V2E, E2E;
        VectorXb boundary, nonManifold;
        if (stats.mMaximumEdgeLength*2 > scale || stats.mMaximumEdgeLength > stats.mAverageEdgeLength * 2) {
            cout << "Input mesh is too coarse for the desired output edge length "
                    "(max input mesh edge length=" << stats.mMaximumEdgeLength
                 << "), subdividing .." << endl;
            build_dedge(F, V, V2E, E2E, boundary, nonManifold);
            subdivide(F, V, V2E, E2E, Feat, boundary, nonManifold, std::min(scale/2, (Float) stats.mAverageEdgeLength*2), deterministic);
            cout << V.cols() << endl;
        }

        /* Compute a directed edge data structure */
        build_dedge(F, V, V2E, E2E, boundary, nonManifold);

        /* Compute adjacency matrix */
        adj = generate_adjacency_matrix_uniform(F, V2E, E2E, nonManifold);

        /* Compute vertex/crease normals */
        if (creaseAngle >= 0)
            generate_crease_normals(F, V, V2E, E2E, boundary, nonManifold, creaseAngle, N, crease_in);
        else
            generate_smooth_normals(F, V, V2E, E2E, nonManifold, N);

        /* Compute dual vertex areas */
        compute_dual_vertex_areas(F, V, V2E, E2E, nonManifold, A);

        mRes.setE2E(std::move(E2E));
    }

    /* Build multi-resolution hierarrchy */
    mRes.setAdj(std::move(adj));
    mRes.setF(std::move(F));
    mRes.setV(std::move(V));
    mRes.setA(std::move(A));
    mRes.setN(std::move(N));
    mRes.setScale(scale);
    mRes.build(deterministic);
    mRes.resetSolution();

    if (align_to_boundaries && !pointcloud) {
        mRes.clearConstraints();
        for (uint32_t i=0; i<3*mRes.F().cols(); ++i) {
            if (mRes.E2E()[i] == INVALID) {
                uint32_t i0 = mRes.F()(i%3, i/3);
                uint32_t i1 = mRes.F()((i+1)%3, i/3);
                Vector3f p0 = mRes.V().col(i0), p1 = mRes.V().col(i1);
                Vector3f edge = p1-p0;
                if (edge.squaredNorm() > 0) {
                    edge.normalize();
                    mRes.CO().col(i0) = p0;
                    mRes.CO().col(i1) = p1;
                    mRes.CQ().col(i0) = mRes.CQ().col(i1) = edge;
                    mRes.CQw()[i0] = mRes.CQw()[i1] = mRes.COw()[i0] =
                        mRes.COw()[i1] = 1.0f;
                }
            }
        }

        for (uint32_t i = 0; i < mRes.F().cols(); ++i) {
            for (int j = 0; j < 3; ++j) {
                if (Feat(j, i) == 1) {
                    uint32_t i0 = mRes.F()(j, i);
                    uint32_t i1 = mRes.F()((j+1)%3, i);
                    Vector3f p0 = mRes.V().col(i0), p1 = mRes.V().col(i1);
                    Vector3f edge = p1-p0;
                    if (mRes.CQw()[i0] > 0 || mRes.CQw()[i1] > 0 ||
                        mRes.COw()[i0] > 0 || mRes.COw()[i1] > 0) {
                        continue;
                    }
                    if (edge.squaredNorm() > 0) {
                        edge.normalize();
                        mRes.CO().col(i0) = p0;
                        mRes.CO().col(i1) = p1;
                        mRes.CQ().col(i0) = mRes.CQ().col(i1) = edge;
                        mRes.CQw()[i0] = mRes.CQw()[i1] = mRes.COw()[i0] =
                            mRes.COw()[i1] = 1.0f;
                    }
                }
            }
        }
        
        mRes.propagateConstraints(rosy, posy);
    }

    if (bvh) {
        bvh->setData(&mRes.F(), &mRes.V(), &mRes.N());
    } else if (smooth_iter > 0) {
        bvh = new BVH(&mRes.F(), &mRes.V(), &mRes.N(), stats.mAABB);
        bvh->build();
    }

    cout << "Preprocessing is done. (total time excluding file I/O: "
         << timeString(timer.reset()) << ")" << endl;

    Optimizer optimizer(mRes, false);
    optimizer.setRoSy(rosy);
    optimizer.setPoSy(posy);
    optimizer.setExtrinsic(extrinsic);

    cout << "Optimizing orientation field .. ";
    cout.flush();
    optimizer.optimizeOrientations(-1);
    optimizer.notify();
    optimizer.wait();
    cout << "done. (took " << timeString(timer.reset()) << ")" << endl;

    std::map<uint32_t, uint32_t> sing;
    compute_orientation_singularities(mRes, sing, extrinsic, rosy);
    cout << "Orientation field has " << sing.size() << " singularities." << endl;
    timer.reset();

    cout << "Optimizing position field .. ";
    cout.flush();
    optimizer.optimizePositions(-1);
    optimizer.notify();
    optimizer.wait();
    cout << "done. (took " << timeString(timer.reset()) << ")" << endl;
    
    //std::map<uint32_t, Vector2i> pos_sing;
    //compute_position_singularities(mRes, sing, pos_sing, extrinsic, rosy, posy);
    //cout << "Position field has " << pos_sing.size() << " singularities." << endl;
    //timer.reset();

    optimizer.shutdown();

    MatrixXf O_extr, N_extr, Nf_extr;
    std::vector<std::vector<TaggedLink>> adj_extr;
    extract_graph(mRes, extrinsic, rosy, posy, adj_extr, O_extr, N_extr,
                  crease_in, crease_out, deterministic);

    MatrixXu F_extr;
    extract_faces(adj_extr, O_extr, N_extr, Nf_extr, F_extr, posy,
            mRes.scale(), crease_out, true, pure_quad, bvh, smooth_iter);
    cout << "Extraction is done. (total time: " << timeString(timer.reset()) << ")" << endl;

    if (output.size() > 0) {
        write_mesh(output, F_extr, O_extr, MatrixXf(), Nf_extr);
    } else {
        verts.resize(O_extr.cols());
        for (int i = 0; i < O_extr.cols(); ++i) {
            verts[i] = std::vector<float>(3);
            for (int j = 0; j < 3; ++j)
                verts[i][j] = O_extr(j, i);
        }

        faces.resize(F_extr.cols());
        std::map<uint32_t, std::pair<uint32_t, std::map<uint32_t, uint32_t>>> irregular;
        size_t nIrregular = 0;
        for (int f = 0; f < F_extr.cols(); ++f) {
            faces[f] = std::vector<int>{};
            if (F_extr.rows() == 4) {
                if (F_extr(2, f) == F_extr(3, f)) {
                    nIrregular++;
                    auto &value = irregular[F_extr(2, f)];
                    value.first = f;
                    value.second[F_extr(0, f)] = F_extr(1, f);
                    continue;
                }
            }
            for (int j = 0; j < F_extr.rows(); ++j)
                faces[f].push_back(F_extr(j, f));
        }
        for (auto item : irregular) {
            auto face = item.second;
            uint32_t v = face.second.begin()->first, first = v, i = 0;
            while (true) {
                faces[face.first].push_back(v);
                v = face.second[v];
                if (v == first || ++i == face.second.size())
                    break;
            }
        }
        cout << "V=" << O_extr.cols() << ", F=" << F_extr.cols();
        if (irregular.size() > 0)
            cout << "(" << irregular.size() << " irregular faces)";
        cout << "." << endl;
        cout.flush();
    }
    if (bvh)
        delete bvh;
}
