// All mesh lists are in ascending cell ID order. Indices refer to that snapshot.
#pragma once

#include "mesh.h"

namespace vertex3d {
using MeshList = std::vector<const Mesh*>;
using CellPair = std::pair<int,int>;
std::vector<CellPair> contact_pairs(const MeshList&, double range, double bin_width);
std::vector<CellPair> exhaustive_pairs(const MeshList&, double range);

struct Contribution { int vertex; Vec3 force; };
struct PairForces {
    std::vector<Contribution> first, second;
    std::size_t contacts = 0;
};
// The default magnitude per unit area is pressure*(1-gap/range), which is
// bounded, so a compressed cell can push through a neighbor. A barrier law
// pressure*(range/gap-1) grows without bound as the gap closes.
PairForces contact_forces(const Mesh&, const Mesh&, double range, double pressure,
                          bool box_filter = true, bool barrier = false);
double contact_law(double gap, double range, bool barrier);
bool surfaces_intersect(const Mesh&, const Mesh&);
bool bodies_intersect(const Mesh&, const Mesh&);

// old and trial meshes have identical connectivity and linearly moving vertices.
// The list overload treats the last argument as an enclosing, stationary wall.
bool swept_crossing(const Mesh& a0, const Mesh& a1, const Mesh& b0, const Mesh& b1);
bool swept_crossing(const MeshList& old_cells, const MeshList& trial_cells, const Mesh& wall);
}
