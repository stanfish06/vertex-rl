#include "solver.h"
#include "contact.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <set>
#include <stdexcept>

namespace vertex3d {
namespace {
template<class F>
void parallel_tasks(int count, int threads, F&& task) {
    std::vector<std::exception_ptr> errors(count);
    #pragma omp parallel for num_threads(threads) schedule(static) if(threads > 1 && count > 1)
    for (int i = 0; i < count; ++i) {
        try { task(i); }
        catch (...) { errors[i] = std::current_exception(); }
    }
    for (auto& error : errors) if (error) std::rethrow_exception(error);
}

MeshList meshes(const std::vector<Cell>& cells) {
    MeshList result;
    for (const auto& cell : cells) result.push_back(&cell.mesh);
    return result;
}

double uniform(std::mt19937& rng) {
    return std::generate_canonical<double, 53>(rng);
}
}

Solver::Solver(Experiment experiment, unsigned seed, int threads)
    : wall_(std::move(experiment.wall)), regions_(std::move(experiment.regions)),
      cells_(std::move(experiment.cells)), material_(experiment.material), random_(seed), threads_(threads) {
    if (threads < 1) throw std::invalid_argument("threads must be positive");
    const auto& m = material_;
    for (double value : {m.volume_stiffness, m.surface_tension, m.pressure, m.growth_half_size, m.smoothing}) {
        if (!std::isfinite(value) || value < 0) throw std::invalid_argument("invalid mechanical coefficient");
    }
    for (double value : {m.contact_range, m.bin_width, m.max_displacement, m.max_step, m.min_step}) {
        if (!std::isfinite(value) || value <= 0) throw std::invalid_argument("invalid length or step bound");
    }
    if (m.min_step > m.max_step) throw std::invalid_argument("min_step exceeds max_step");
    if (m.remesh_level < 0 || m.remesh_level > 6) throw std::invalid_argument("remesh_level must be 0 to 6");
    validate(wall_);
    if (regions_.empty()) regions_.push_back({Zone{0, [](Vec3) { return true; }}});
    std::set<int> ids;
    for (const Region& r : regions_) {
        if (!ids.insert(r.zone.id).second) throw std::invalid_argument("duplicate zone ID");
        if (!std::isfinite(r.growth_rate) || r.growth_rate < 0 ||
            !std::isfinite(r.division_rate) || r.division_rate < 0 || !(r.division_volume > 0) ||
            !(r.target_size > 0) || !std::isfinite(r.death_rate) || r.death_rate < 0) {
            throw std::invalid_argument("invalid zone growth or division settings");
        }
        zones_.push_back(r.zone);
    }
    for (Cell& cell : cells_) {
        validate(cell.mesh);
        cell.id = next_id_++;
        measure_cell(cell);
        if (cell.target_volume == 0) cell.target_volume = cell.geometry.volume;
        if (!(cell.target_volume > 0) || !std::isfinite(cell.target_volume)) {
            throw std::invalid_argument("target_volume must be positive and finite");
        }
        if (!contains(wall_, cell.mesh)) throw std::invalid_argument("cell extends outside the wall");
    }
    auto list = meshes(cells_);
    for (auto [i, j] : contact_pairs(list, m.contact_range, m.bin_width)) {
        if (bodies_intersect(*list[i], *list[j])) throw std::invalid_argument("initial cell meshes intersect");
    }
}

void Solver::measure_cell(Cell& cell) const {
    cell.geometry = measure(cell.mesh);
    if (!(cell.geometry.volume > 0) || !std::isfinite(cell.geometry.volume) ||
        !std::isfinite(cell.geometry.area) || !std::isfinite(norm(cell.geometry.centroid))) {
        throw std::runtime_error("invalid cell volume, area or centroid");
    }
    cell.zone = classify(cell.geometry.centroid, zones_);
}

const Region& Solver::region(int id) const {
    return *std::find_if(regions_.begin(), regions_.end(), [id](const Region& r) { return r.zone.id == id; });
}

bool Solver::rebuild(Cell& cell) const {
    // The rebuilt vertices lie on the current surface along rays from the
    // centroid. Chords of a dented surface can bulge past it, so contract toward
    // the centroid until every new face centroid lies inside the current cell.
    // The rebuilt cell then occupies only space the cell already occupied.
    Mesh sphere_mesh;
    try { sphere_mesh = star_remesh(cell.mesh, cell.geometry.centroid, material_.remesh_level); }
    catch (const std::invalid_argument&) { return false; }
    for (double scale : {0.99, 0.97, 0.94, 0.9}) {
        Cell copy = cell;
        copy.mesh = sphere_mesh;
        for (Vec3& v : copy.mesh.vertices) v = cell.geometry.centroid + (v - cell.geometry.centroid) * scale;
        bool inside = true;
        for (const Face& f : copy.mesh.faces) {
            Vec3 c = (copy.mesh.vertices[f[0]] + copy.mesh.vertices[f[1]] + copy.mesh.vertices[f[2]]) / 3;
            if (locate_point(cell.mesh, c) < 0) { inside = false; break; }
        }
        if (!inside || !contains(wall_, copy.mesh)) continue;
        try { measure_cell(copy); } catch (const std::exception&) { return false; }
        cell = std::move(copy);
        return true;
    }
    return false;
}

bool Solver::cut(std::vector<Cell>& tissue, int id, Vec3 normal, int& next_id, double born) const {
    auto parent = std::find_if(tissue.begin(), tissue.end(), [id](const Cell& c) { return c.id == id; });
    if (parent == tissue.end()) throw std::invalid_argument("unknown cell ID " + std::to_string(id));
    if (!region(parent->zone).allow_division) return false;
    Cell a, b;
    try {
        auto halves = split(parent->mesh, parent->geometry.centroid, normal);
        a.mesh = std::move(halves.first);
        b.mesh = std::move(halves.second);
        measure_cell(a);
        measure_cell(b);
    } catch (const std::invalid_argument&) {
        return false;
    }
    double total = a.geometry.volume + b.geometry.volume;
    if (std::abs(total - parent->geometry.volume) > 1e-9 * parent->geometry.volume) return false;
    if (material_.remesh_level > 0) {
        // Daughters inherit the cut triangulation, including the long thin
        // triangles of the planar cap. Rebuild each on a uniform sphere.
        rebuild(a);
        rebuild(b);
        total = a.geometry.volume + b.geometry.volume;
    }
    a.born = b.born = born;
    a.growing = b.growing = parent->growing;
    a.target_volume = parent->target_volume * (a.geometry.volume / total);
    b.target_volume = parent->target_volume - a.target_volume;
    if (!(a.target_volume > 0 && b.target_volume > 0)) return false;
    a.id = next_id++;
    b.id = next_id++;
    tissue.erase(parent);
    tissue.push_back(std::move(a));
    tissue.push_back(std::move(b));
    return true;
}

bool Solver::divide(int id, Vec3 normal) {
    if (!std::isfinite(norm(normal)) || !(norm(normal) > 0)) throw std::invalid_argument("invalid plane normal");
    auto proposed = cells_;
    int next = next_id_;
    if (!cut(proposed, id, normal, next, time_)) return false;
    cells_.swap(proposed);
    next_id_ = next;
    return true;
}

std::vector<std::vector<Vec3>> Solver::forces() const {
    const Material& m = material_;
    auto list = meshes(cells_);
    auto pairs = contact_pairs(list, m.contact_range, m.bin_width);
    std::vector<std::vector<Vec3>> result(cells_.size());
    parallel_tasks(cells_.size(), threads_, [&](int i) {
        result[i] = internal_forces(*list[i], cells_[i].target_volume, m.volume_stiffness, m.surface_tension);
        // Without a wall reaction a cell growing into the wall can only shorten
        // its step until the minimum and fail.
        auto wall = wall_forces(*list[i], wall_, m.contact_range, m.pressure, m.contact_barrier);
        for (size_t k = 0; k < wall.size(); ++k) result[i][k] += wall[k];
        if (m.smoothing > 0) {
            // Umbrella Laplacian: keeps vertices evenly spread and damps the
            // short-wavelength folds that contact pressure excites at corners.
            const Mesh& mesh = *list[i];
            std::vector<Vec3> sum(mesh.vertices.size());
            std::vector<int> degree(mesh.vertices.size(), 0);
            for (const auto& item : mesh.edges) {
                auto [p, q] = item.first;
                sum[p] += mesh.vertices[q]; sum[q] += mesh.vertices[p];
                ++degree[p]; ++degree[q];
            }
            for (size_t k = 0; k < sum.size(); ++k)
                if (degree[k]) result[i][k] += (sum[k] / degree[k] - mesh.vertices[k]) * m.smoothing;
        }
    });
    std::vector<PairForces> buffers(pairs.size());
    parallel_tasks(pairs.size(), threads_, [&](int k) {
        auto [a,b] = pairs[k];
        buffers[k] = contact_forces(*list[a], *list[b], m.contact_range, m.pressure, true, m.contact_barrier);
    });
    // Cells stay in ascending ID order. The grid supplies sorted unique pairs.
    std::vector<std::vector<int>> incident(cells_.size());
    for (size_t k = 0; k < pairs.size(); ++k) {
        incident[pairs[k].first].push_back(k);
        incident[pairs[k].second].push_back(k);
    }
    parallel_tasks(cells_.size(), threads_, [&](int i) {
        for (int k : incident[i]) {
            const auto& contributions = pairs[k].first == i ? buffers[k].first : buffers[k].second;
            for (const auto& term : contributions) result[i][term.vertex] += term.force;
        }
    });
    return result;
}

void Solver::configure(std::vector<Region> regions) {
    // The driver validates parameter assignments before replacing the list.
    regions_.swap(regions);
}

void Solver::policy_step(double deadline, std::vector<Action> actions) {
    std::sort(actions.begin(), actions.end(), [](const Action& a, const Action& b) { return a.id < b.id; });
    if (actions.size() != cells_.size()) throw std::invalid_argument("expected exactly one action per live cell");
    for (size_t i = 0; i < cells_.size(); ++i) {
        if (cells_[i].id != actions[i].id) throw std::invalid_argument("missing, duplicate or unknown action ID");
    }
    auto proposed = cells_;
    int next = next_id_;
    for (const auto& action : actions) {
        auto cell = std::lower_bound(proposed.begin(), proposed.end(), action.id,
            [](const Cell& c, int id) { return c.id < id; });
        const Region& r = region(cell->zone);
        switch (action.fate) {
        case Fate::grow: if (r.allow_growth) cell->growing = true; break;
        case Fate::stop: if (r.allow_stop) cell->growing = false; break;
        case Fate::divide: cut(proposed, cell->id, action.normal, next, time_); break;
        case Fate::die: if (r.allow_death) proposed.erase(cell); break;
        }
    }
    double start = time_;
    std::swap(next_id_, next);
    cells_.swap(proposed);
    try { advance(deadline, false); }
    catch (...) {
        if (time_ == start) {
            cells_.swap(proposed);
            next_id_ = next;
        }
        throw;
    }
}

bool Solver::repair(const std::vector<Cell>& trial) {
    std::vector<bool> offending(cells_.size(), false);
    auto old = meshes(cells_);
    for (size_t i = 0; i < cells_.size(); ++i) {
        try { validate_geometry(trial[i].mesh); } catch (const std::invalid_argument&) { offending[i] = true; }
        if (swept_crossing(cells_[i].mesh, trial[i].mesh, wall_, wall_)) offending[i] = true;
    }
    for (auto [i, j] : contact_pairs(old, material_.contact_range + material_.max_displacement, material_.bin_width)) {
        if (swept_crossing(cells_[i].mesh, trial[i].mesh, cells_[j].mesh, trial[j].mesh)) offending[i] = offending[j] = true;
    }
    bool repaired = false;
    for (size_t i = 0; i < cells_.size(); ++i) if (offending[i]) repaired |= rebuild(cells_[i]);
    return repaired;
}

void Solver::advance(double deadline, bool builtin_fates) {
    if (!std::isfinite(deadline) || deadline <= time_) throw std::invalid_argument("deadline must be later than time");
    int repairs = 0;
    while (time_ < deadline && !cells_.empty()) {
        bool restart = false;
        auto force = forces();
        double max_speed = 0;
        for (const auto& cell : force) for (Vec3 f : cell) {
            double speed = norm(f);
            if (!std::isfinite(speed)) throw std::runtime_error("nonfinite vertex force");
            max_speed = std::max(max_speed, speed);
        }
        double remaining = deadline - time_;
        double dt = std::min(remaining, material_.max_step);
        if (max_speed > 0) dt = std::min(dt, material_.max_displacement / max_speed);
        auto proposed = cells_;
        auto old = meshes(cells_);
        auto trial = meshes(proposed);
        std::uint64_t rejected = 0;
        while (true) {
            if (!(dt > 0) || time_ + dt == time_) throw std::runtime_error("step cannot advance the clock");
            for (size_t i = 0; i < cells_.size(); ++i) {
                for (size_t j = 0; j < force[i].size(); ++j) {
                    proposed[i].mesh.vertices[j] = cells_[i].mesh.vertices[j] + force[i][j]*dt;
                }
            }
            bool valid = !swept_crossing(old, trial, wall_);
            if (valid) {
                try {
                    for (auto& cell : proposed) {
                        // Connectivity is immutable during motion. A cell can still
                        // collapse, invert or self-intersect on an excessive trial.
                        validate_geometry(cell.mesh);
                        measure_cell(cell);
                    }
                } catch (const std::invalid_argument&) { valid = false; }
            }
            if (valid) break;
            if (dt <= material_.min_step) {
                // Repair the cells whose trial folds or crosses by rebuilding
                // their meshes, then recompute forces. Each repair contracts a
                // cell slightly, so repeated repairs cannot cycle indefinitely.
                if (material_.remesh_level > 0 && repairs < 8 && repair(proposed)) { ++repairs; restart = true; break; }
                throw std::runtime_error("crossing or invalid mesh at minimum step");
            }
            dt = std::max(material_.min_step, dt*0.5);
            ++rejected;
        }

        if (restart) continue;
        // Fate sampling follows accepted geometry in stable cell ID order.
        // The working generator and ID counter commit with the cells and time.
        auto rng = random_;
        int next = next_id_;
        std::vector<std::pair<int, Vec3>> divisions;
        std::vector<int> deaths;
        for (auto& cell : proposed) {
            const Region& r = region(cell.zone);
            if (cell.growing && r.allow_growth && cell.target_volume < r.target_size) {
                double rate = r.growth_rate;
                if (material_.growth_half_size > 0) rate /= 1 + cell.geometry.volume / material_.growth_half_size;
                cell.target_volume = std::min(r.target_size, cell.target_volume + rate * dt);
            }
            if (!std::isfinite(cell.target_volume)) throw std::runtime_error("growth produced a nonfinite target");
            if (!builtin_fates) continue;
            if (r.allow_death && r.death_rate > 0 && uniform(rng) < -std::expm1(-r.death_rate*dt)) {
                deaths.push_back(cell.id);
                continue;
            }
            if (r.allow_division && r.division_rate > 0 && cell.geometry.volume >= r.division_volume &&
                uniform(rng) < -std::expm1(-r.division_rate*dt)) {
                double z = 2*uniform(rng)-1, phi = 2*std::acos(-1.)*uniform(rng);
                double radius = std::sqrt(std::max(0., 1-z*z));
                divisions.emplace_back(cell.id, Vec3{radius*std::cos(phi), radius*std::sin(phi), z});
            }
        }
        if (!deaths.empty()) proposed.erase(std::remove_if(proposed.begin(), proposed.end(), [&](const Cell& c) {
            return std::binary_search(deaths.begin(), deaths.end(), c.id);
        }), proposed.end());
        for (const auto& d : divisions) cut(proposed, d.first, d.second, next, time_ + dt);
        cells_.swap(proposed);
        random_ = rng;
        next_id_ = next;
        time_ = dt == remaining ? deadline : time_ + dt;
        retries_ += rejected;
        ++accepted_;
    }
}
}
