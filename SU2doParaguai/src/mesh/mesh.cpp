#include "mesh/mesh.hpp"

#include <iostream>
#include <utility>
#include <map>
#include <algorithm>
#include <limits>
#include <vector>
#include <cmath>

namespace FVMesh
{
    CMesh::CMesh()
        : nDime(0), nElem(0), nPoin(0), nMark(0)
    {
    }

    CMesh::CMesh(int dime, unsigned long elem, unsigned long poin, unsigned long mark, std::vector<FVUtil::Vec2>&& nodes, std::vector<Elem>&& elems, std::vector<Marker>&& marks, bool isStatic_)
        : Nodes(std::move(nodes)), Elems(std::move(elems)), Marks(std::move(marks)), nDime(dime), nElem(elem), nPoin(poin), nMark(mark)
    {
        isStatic = isStatic_;
        Centroids.reserve(nElem);

        for (size_t i = 0; i < nPoin; i++)
            Volumes.push_back(0.0);

        for (size_t i = 0; i < Elems.size(); i++)
        {
            if (Elems[i].type == ElemType::TRIANGLE)
            {
                const Elem& e = Elems[i];
                const FVUtil::Vec2 centr =
                    (Nodes[e.nodes[0]] + Nodes[e.nodes[1]] +
                        Nodes[e.nodes[2]]) * (1.0 / 3.0);
                Centroids.push_back(centr);
            }
            else
            {
                const Elem& e = Elems[i];
                double area = 0.0, x_sum = 0.0, y_sum = 0.0;

                for (size_t j = 0; j < e.nodes.size(); j++)
                {
                    const FVUtil::Vec2 pos_i = Nodes[e.nodes[j]];
                    const FVUtil::Vec2 pos_nex =
                        Nodes[e.nodes[(j + 1) % e.nodes.size()]];
                    const double cross =
                        pos_i.x * pos_nex.y - pos_nex.x * pos_i.y;

                    area += cross;
                    x_sum += (pos_i.x + pos_nex.x) * cross;
                    y_sum += (pos_i.y + pos_nex.y) * cross;
                }

                area *= 0.5;
                if (std::abs(area) <= std::numeric_limits<double>::epsilon())
                {
                    std::cerr << "ERROR: Degenerate element " << i
                        << " while calculating its centroid.\n";
                    Centroids.emplace_back(0.0, 0.0);
                }
                else
                {
                    const double factor = 1.0 / (6.0 * area);
                    Centroids.emplace_back(x_sum * factor, y_sum * factor);
                }
            }
        }
    }

    void CMesh::InitDualMesh()
    {
        Edges.clear();
        BoundaryEdges.clear();
        nIntEdges = 0;
        nBoundEdges = 0;
        std::fill(Volumes.begin(), Volumes.end(), 0.0);

        std::map<std::pair<int, int>, size_t> edge_map;

        // Build primal edges and their adjacent-element connectivity.
        for (size_t i = 0; i < Elems.size(); i++)
        {
            const auto& nodes = Elems[i].nodes;
            const size_t num_nodes = nodes.size();

            for (size_t k = 0; k < num_nodes; k++)
            {
                const int n1 = nodes[k];
                const int n2 = nodes[(k + 1) % num_nodes];
                const std::pair<int, int> edge_key = {
                    std::min(n1, n2), std::max(n1, n2)
                };

                const auto iter = edge_map.find(edge_key);
                if (iter == edge_map.end())
                {
                    Edge new_edge;
                    new_edge.nBeg = n1;
                    new_edge.nEnd = n2;
                    new_edge.eLef = static_cast<int>(i);
                    new_edge.eRig = -1;
                    Edges.push_back(new_edge);
                    edge_map[edge_key] = Edges.size() - 1;
                }
                else
                {
                    Edge& shared_edge = Edges[iter->second];
                    if (shared_edge.eRig != -1)
                    {
                        std::cerr << "ERROR: Non-manifold edge shared by "
                            "more than two elements: "
                            << n1 << " - " << n2 << '\n';
                    }
                    else
                    {
                        shared_edge.eRig = static_cast<int>(i);
                    }
                }
            }
        }

        for (auto& e : Edges)
            e.IsBoundary = (e.eRig == -1);

        // Calculate integrated median-dual normals and dual volumes.
        for (size_t i = 0; i < Edges.size(); i++)
        {
            Edge& currEdge = Edges[i];
            const int owner = currEdge.nBeg;
            const int neighb = currEdge.nEnd;
            const int leftEl = currEdge.eLef;
            const int rightEl = currEdge.eRig;

            const FVUtil::Vec2 edge_vec = Nodes[neighb] - Nodes[owner];
            const FVUtil::Vec2 midpoint =
                (Nodes[owner] + Nodes[neighb]) * 0.5;

            if (rightEl != -1)
            {
                const FVUtil::Vec2 centroid_vec =
                    Centroids[rightEl] - Centroids[leftEl];

                FVUtil::Vec2 normal(-centroid_vec.y, centroid_vec.x);

                if (FVUtil::Vec2::Dot(normal, edge_vec) < 0.0)
                    normal = normal * (-1.0);

                currEdge.Normal = normal;
            }
            else
            {
                // FIX: Calculate the internal dual face for boundary edges 
                // connecting the element centroid to the edge midpoint.
                const FVUtil::Vec2 face_vec = midpoint - Centroids[leftEl];
                FVUtil::Vec2 normal(-face_vec.y, face_vec.x);

                if (FVUtil::Vec2::Dot(normal, edge_vec) < 0.0)
                    normal = normal * (-1.0);

                currEdge.Normal = normal;
            }

            Volumes[owner] += FVUtil::Vec2::TriArea(
                Centroids[leftEl] - Nodes[owner],
                midpoint - Nodes[owner]);

            Volumes[neighb] += FVUtil::Vec2::TriArea(
                Centroids[leftEl] - Nodes[neighb],
                midpoint - Nodes[neighb]);

            if (rightEl != -1)
            {
                Volumes[owner] += FVUtil::Vec2::TriArea(
                    Centroids[rightEl] - Nodes[owner],
                    midpoint - Nodes[owner]);

                Volumes[neighb] += FVUtil::Vec2::TriArea(
                    Centroids[rightEl] - Nodes[neighb],
                    midpoint - Nodes[neighb]);
            }
        }

        auto boundary_start_iter = std::partition(
            Edges.begin(), Edges.end(),
            [](const Edge& e) { return !e.IsBoundary; });

        nIntEdges = static_cast<ulong>(
            std::distance(Edges.begin(), boundary_start_iter));
        nBoundEdges = static_cast<ulong>(Edges.size()) - nIntEdges;
        BoundaryEdges.reserve(nBoundEdges);

        for (ulong e = nIntEdges; e < Edges.size(); e++)
        {
            const Edge& currEdge = Edges[e];
            BoundaryEdge edge;
            edge.nBeg = currEdge.nBeg;
            edge.nEnd = currEdge.nEnd;

            const FVUtil::Vec2 edge_vec =
                Nodes[edge.nEnd] - Nodes[edge.nBeg];
            const FVUtil::Vec2 midpoint =
                (Nodes[edge.nBeg] + Nodes[edge.nEnd]) * 0.5;

            // Start with the right-hand normal
            FVUtil::Vec2 normal(edge_vec.y, -edge_vec.x);

            // Find the edge's orientation in the element's local node list
            const auto& nodes = Elems[currEdge.eLef].nodes;
            bool is_forward = false;
            for (size_t k = 0; k < nodes.size(); k++)
            {
                if (nodes[k] == edge.nBeg && nodes[(k + 1) % nodes.size()] == edge.nEnd)
                {
                    is_forward = true;
                    break;
                }
            }

            // In a CCW element, if the edge connects nBeg -> nEnd sequentially, 
            // the interior is to the left, making the right-hand normal outward.
            // If the sequence is reversed, we must flip the normal.
            if (!is_forward)
                normal = normal * (-1.0);

            edge.Normal = normal;
            BoundaryEdges.push_back(edge);
        }
        // DEBUG TEST 1: internal normal points from nBeg to nEnd.
        double min_orientation_dot = std::numeric_limits<double>::max();
        ulong bad_orientation_count = 0;

        // FIX: Changed from nIntEdges to Edges.size()
        for (ulong e = 0; e < Edges.size(); e++)
        {
            const Edge& currEdge = Edges[e];
            const FVUtil::Vec2 edge_vec =
                Nodes[currEdge.nEnd] - Nodes[currEdge.nBeg];
            const double orientation_dot =
                FVUtil::Vec2::Dot(currEdge.Normal, edge_vec);

            min_orientation_dot =
                std::min(min_orientation_dot, orientation_dot);
            if (orientation_dot <= 0.0)
                bad_orientation_count++;
        }

        std::cout << "Mesh Initiated\n\tEdges:\n"
            << "\tTotal Edges: " << Edges.size() << "\n"
            << "\tInternal Edges: " << nIntEdges << "\n"
            << "\tBoundary Edges: " << nBoundEdges << "\n";

        std::cout << "\tDEBUG - Internal normal orientation:\n"
            << "\t\tMinimum S_ij dot (x_j - x_i): "
            << min_orientation_dot << "\n"
            << "\t\tInvalid internal normals: "
            << bad_orientation_count << "\n";

        // DEBUG TEST 2: sum of integrated face normals around each dual CV.
        std::vector<FVUtil::Vec2> closure(
            nPoin, FVUtil::Vec2(0.0, 0.0));

        // FIX: Changed from nIntEdges to Edges.size()
        for (ulong e = 0; e < Edges.size(); e++)
        {
            const Edge& currEdge = Edges[e];
            closure[currEdge.nBeg] =
                closure[currEdge.nBeg] + currEdge.Normal;
            closure[currEdge.nEnd] =
                closure[currEdge.nEnd] + currEdge.Normal * (-1.0);
        }

        for (ulong e = 0; e < nBoundEdges; e++)
        {
            const BoundaryEdge& currEdge = BoundaryEdges[e];
            const FVUtil::Vec2 half_normal = currEdge.Normal * 0.5;

            closure[currEdge.nBeg] =
                closure[currEdge.nBeg] + half_normal;
            closure[currEdge.nEnd] =
                closure[currEdge.nEnd] + half_normal;
        }

        double max_closure_residual = 0.0;
        double total_normal_magnitude = 0.0;
        ulong normal_count = 0;

        for (ulong i = 0; i < nPoin; i++)
            max_closure_residual = std::max(
                max_closure_residual, FVUtil::Vec2::Mag(closure[i]));

        // FIX: Changed from nIntEdges to Edges.size()
        for (ulong e = 0; e < Edges.size(); e++)
        {
            total_normal_magnitude += FVUtil::Vec2::Mag(Edges[e].Normal);
            normal_count++;
        }

        for (ulong e = 0; e < nBoundEdges; e++)
        {
            total_normal_magnitude +=
                FVUtil::Vec2::Mag(BoundaryEdges[e].Normal);
            normal_count++;
        }

        const double mean_normal_magnitude = normal_count > 0
            ? total_normal_magnitude / static_cast<double>(normal_count)
            : 0.0;

        const double closure_tolerance =
            1.0e-10 * std::max(1.0, mean_normal_magnitude);

        ulong bad_closure_count = 0;
        for (ulong i = 0; i < nPoin; i++)
        {
            if (FVUtil::Vec2::Mag(closure[i]) > closure_tolerance)
                bad_closure_count++;
        }

        std::cout << "\tDEBUG - Dual control-volume normal closure:\n"
            << "\t\tMaximum |sum(S_f)|: "
            << max_closure_residual << "\n"
            << "\t\tMean face-normal magnitude: "
            << mean_normal_magnitude << "\n"
            << "\t\tClosure tolerance: " << closure_tolerance << "\n"
            << "\t\tPoints failing closure: "
            << bad_closure_count << " / " << nPoin << "\n";

        double sum_v = 0.0;
        std::cout << "Cell volumes:\n";

        for (size_t i = 0; i < Nodes.size(); i++)
        {
            sum_v += Volumes[i];
            std::cout << "\tNode (" << Nodes[i].x << ", " << Nodes[i].y
                << "): V = " << Volumes[i] << " | counter: \n";
        }

        std::cout << "Sum of V = " << sum_v << "\n";

        for (size_t i = 0; i < Elems.size(); i++)
            std::cout << "Elem centroid: " << Centroids[i].x << ", "
            << Centroids[i].y << "\n";

        std::cout << "Edges:\n";
        for (size_t i = 0; i < Edges.size(); i++)
        {
            std::cout << "\t" << Edges[i].nBeg << " -> " << Edges[i].nEnd
                << "\n\t\tBoundary: " << Edges[i].IsBoundary
                << " | Lef: " << Edges[i].eLef
                << " | Rig: " << Edges[i].eRig
                << " | Norm: " << Edges[i].Normal.x
                << ", " << Edges[i].Normal.y << "\n";
        }
    }
}
