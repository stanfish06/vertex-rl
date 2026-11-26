#include "solver_2d.h"
#include <iomanip>
#include <sstream>

// One solver lives in each worker process. All C ABI calls run on its request thread.
namespace {
Solver2D solver(true);
std::string last_error;
std::string payload;

std::string state_json() {
    const auto& tissue = solver.observation();
    std::ostringstream out;
    out << std::setprecision(17) << "{\"time\":" << tissue.time << ",\"cells\":[";
    std::array<int, 5> zones{};
    int live = 0;
    double area = 0;
    bool first = true;
    for (const auto& cell : tissue.cells) {
        if (!first) out << ',';
        first = false;
        if (!cell.dead) {
            ++live;
            ++zones[cell.zone];
            area += cell.area;
        }
        out << "{\"id\":" << cell.id << ",\"zone\":" << cell.zone
            << ",\"dead\":" << (cell.dead ? "true" : "false")
            << ",\"growing\":" << (cell.growing ? "true" : "false")
            << ",\"center\":[" << cell.x << ',' << cell.y << "],\"area\":" << cell.area
            << ",\"target_area\":" << cell.target_area << ",\"vertices\":[";
        for (size_t j = 0; j < cell.vertices.size(); ++j) {
            if (j) out << ',';
            out << '[' << cell.vertices[j].x << ',' << cell.vertices[j].y << ']';
        }
        out << "]}";
    }
    out << "],\"measurements\":{\"cell_count\":" << live << ",\"total_area\":" << area
        << ",\"zone_counts\":[";
    for (size_t z = 0; z < zones.size(); ++z) {
        if (z) out << ',';
        out << zones[z];
    }
    out << "]}}";
    return out.str();
}
}

extern "C" {
const char* vertex2d_error() { return last_error.c_str(); }

int vertex2d_reset(unsigned seed, int n, const char** names, const double* values) {
    try {
        // Validate all overrides before replacing an existing episode.
        for (int i = 0; i < n; ++i) {
            auto* option = find_parameter(names[i]);
            float value = static_cast<float>(values[i]);
            bool positive = option && (option->core || option->name == std::string("maxDR") ||
                                        option->name == std::string("lthUP"));
            if (!option || !std::isfinite(value) || value < 0 || (positive && value == 0)) {
                last_error = "invalid model parameter " + std::string(names[i]);
                return -2;
            }
        }
        for (auto& parameter : model_parameters) parameter.reset();
        for (int i = 0; i < n; ++i) find_parameter(names[i])->set(values[i]);
        verbose_solver = false;
        set_solver_threads(1);
        solver.reset(seed);
        return 0;
    } catch (const std::exception& error) {
        last_error = error.what();
        return -1;
    }
}

int vertex2d_advance(double deadline) {
    try {
        solver.advance(deadline);
        return 0;
    } catch (const std::exception& error) {
        last_error = error.what();
        return -1;
    }
}

const char* vertex2d_measure() {
    payload = state_json();
    return payload.c_str();
}

int vertex2d_export(int frame) {
    try {
        export_tissue(solver.observation(), frame);
        return 0;
    } catch (const std::exception& error) {
        last_error = error.what();
        return -1;
    }
}
}
