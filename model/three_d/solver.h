#pragma once
#include "mesh.h"
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

namespace vertex3d {

struct Material {
    double volume_stiffness = 1.0;
    double surface_tension = 0.05;
    double pressure = 1.0;
    double contact_range = 0.1;
    double bin_width = 1.0;
    double max_displacement = 0.05;
    double max_step = 0.01;
    double min_step = 1e-8;
};

struct Region {
    Zone zone;
    bool allow_growth = true;
    bool allow_division = true;
    double growth_rate = 0.0;     // target volume per unit time
    double division_volume = std::numeric_limits<double>::infinity();
    double division_rate = 0.0;   // hazard per unit time after reaching division_volume
};

struct Cell {
    int id = 0;
    Mesh mesh;
    double target_volume = 0;
    bool growing = true;
    int zone = 0;
    Measures geometry{};
};

struct Experiment {
    Mesh wall;
    std::vector<Region> regions;
    std::vector<Cell> cells;
    Material material;
};

class Solver {
public:
    Solver(Experiment experiment, unsigned seed, int threads);
    void advance(double deadline);
    bool divide(int cell_id, Vec3 normal);
    const std::vector<Cell>& cells() const { return cells_; }
    const std::vector<Region>& regions() const { return regions_; }
    double time() const { return time_; }
    std::uint64_t retries() const { return retries_; }
    std::uint64_t accepted_steps() const { return accepted_; }

private:
    Mesh wall_;
    std::vector<Region> regions_;
    std::vector<Zone> zones_;
    std::vector<Cell> cells_;
    Material material_;
    std::mt19937 random_;
    int threads_, next_id_ = 1;
    double time_ = 0;
    std::uint64_t retries_ = 0, accepted_ = 0;
    void measure_cell(Cell& cell) const;
    const Region& region(int id) const;
    bool cut(std::vector<Cell>& tissue, int id, Vec3 normal, int& next_id) const;
    std::vector<std::vector<Vec3>> forces() const;
};
}
