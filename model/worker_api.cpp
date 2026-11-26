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
    std::vector<int> zones(training2d.mode ? training2d.zones.size() : 5);
    int live = 0;
    double area = 0;
    bool first = true;
    for (const auto& cell : tissue.cells) {
        if (training2d.mode && cell.dead) continue;
        if (!first) out << ',';
        first = false;
        if (!cell.dead) {
            ++live;
            ++zones.at(training2d.mode ? training2d.index.at(cell.zone) : cell.zone);
            area += cell.area;
        }
        out << "{\"id\":" << (training2d.mode ? cell.agent_id : cell.id) << ",\"zone\":" << cell.zone;
        if (training2d.mode) out << ",\"age\":" << cell.age;
        out << ",\"dead\":" << (cell.dead ? "true" : "false")
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

struct ZoneSpec {
    int id, rule, flags, positive, count;
    double growth, target, division, death, lower, upper;
    const double* points;
};
int vertex2d_reset_training(unsigned seed, int threads, int n, const char** names, const double* values,
                            int mode, int count, const ZoneSpec* specs) {
    try {
        if (threads < 1) {
            last_error = "thread count must be positive";
            return -2;
        }
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
        Training2D settings;
        settings.mode=mode;
        if (mode) {
            for (int i=0; i<count; ++i) {
                const auto& spec=specs[i];
                Zone2D z;
                z.id=spec.id; z.rule=spec.rule;
                z.allow_growth=spec.flags&1; z.allow_stop=spec.flags&2;
                z.allow_division=spec.flags&4; z.allow_death=spec.flags&8;
                z.growth_rate=spec.growth; z.target_size=spec.target;
                z.division_rate=spec.division; z.death_rate=spec.death;
                z.lower=spec.lower; z.upper=spec.upper; z.positive=spec.positive;
                for (int j=0; j<spec.count; ++j) z.points.push_back({spec.points[2*j],spec.points[2*j+1]});
                if (z.rule==5 && !simple_ring(z.points)) {
                    last_error="invalid closed zone polygon";
                    return -2;
                }
                settings.index.emplace(z.id,settings.zones.size());
                settings.zones.push_back(std::move(z));
            }
        }
        for (auto& parameter : model_parameters) parameter.reset();
        for (int i=0; i<n; ++i) find_parameter(names[i])->set(values[i]);
        for (auto& z : settings.zones) if (!z.rule && z.target_size==0) z.target_size=1.5*CV[z.id];
        training2d=std::move(settings);
        verbose_solver = false;
        set_solver_threads(threads);
        solver.reset(seed);
        return 0;
    } catch (const std::exception& error) {
        last_error = error.what();
        return -1;
    }
}

int vertex2d_reset(unsigned seed, int threads, int n, const char** names, const double* values) {
    return vertex2d_reset_training(seed,threads,n,names,values,0,0,nullptr);
}

void vertex2d_control(int id, int parameter, double value) {
    auto& z=training2d.zones[training2d.index.at(id)];
    if (parameter==0) z.growth_rate=value;
    else if (parameter==1) z.target_size=value;
    else if (parameter==2) z.division_rate=value;
    else z.death_rate=value;
}

int vertex2d_actions(int count, const int* agents, const int* actions, const double* angles) {
    try {
        std::map<int,int> slots;
        for (const auto& cell : solver.observation().cells) if (!cell.dead) slots.emplace(cell.agent_id,cell.id);
        for (int i=0; i<count; ++i) solver.act(slots.at(agents[i]),actions[i],angles[i]);
        return 0;
    } catch (const std::exception& error) {
        last_error=error.what();
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
