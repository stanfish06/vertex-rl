#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

//****************INPUT PARAMETERS************************
//MODEL PARAMETERS
float k0=99;      //cell-cell repulsion
float kS=2;      //ECM shoulder repulsion
float lS=0.5;      //ECM shoulder thickness
float kC=0.5;
//float e_g0=0.25;    //cell line tension (should be below 0.434313)
float Atarget=1;   //cell target area (should be kept 1=normalized)
float kA=1.2; //1.2       //cell-area elastic modulus (should be kept 1=normalized)
float e_g0=0.7; //0.4    //cell line tension (should be below 0.434313)
// Threshold of line tension should be calculated by balancing with area compression, slightly higher than threshold is fine as desired area will change.
// However, if initial area is too small and cell's desired area changes slow, cell area might shrink to 0 as there is no nonzero energy minimum
float tMAX=2000;   //total integration time
float outFreq=0.5;   //frequency of outputing
float h0=0.001;    //time step (should not be too big)
float L=4;        //width of the simulation box

//PARAMETERS FOR EDGE MANIPULATION
float lthUP=0.2;   //if edge longer than this, it gets divided
float lthDOWN=0.01;//if edge shorter than this, its vertices are merged
float lALERT=0.4;  //width of the region around edges, where vertex-edge interactions are computed
float lDEEP=0.5;
//OTHER PARAMETERS
float h;
float Time;
float maxDR=0.001;
float ddTH, ddTHS, vMAX;
float TimeCount=0;
int kount=0;
unsigned int placement_seed = 646;
bool reference_output = false;
bool verbose_solver = true;

// Parallelism parameters
int num_threads_count = 16;
int num_threads_division = 16;
int num_threads_edge_manipulation = 10;
int num_threads_force = 16;
int num_threads_calculation = 16;
int num_threads_update = 16;

float scale_factor = 2.0;

// Zone interfaces
float proliferative_zone_interface[] = {7 * scale_factor, 11 * scale_factor};
float hypertrophic_zone_interface[] = {4 * scale_factor, 7 * scale_factor};
float resting_zone_interface[] = {11 * scale_factor, (14 + 2 * 2 * sqrt(2)) * scale_factor};
float ossification_interface[] = {0 * scale_factor, 4 * scale_factor};
float pericondrium_thickness = 1.0;

float proliferative_zone_area = (proliferative_zone_interface[1] - proliferative_zone_interface[0]) * 4 * scale_factor;
// Polygon interface (all subpolygons have to be convex)
float polygon_interface_y = 10 * scale_factor;

// Region ECM (scale factor for lS)
float proliferative_ECM = 2;
float hypertrophic_ECM = 3;
float pericondrium_ECM = 1;
float ossification_ECM = 3;
float resting_ECM = 1;

float ECM[] = {
    resting_ECM,
    proliferative_ECM,
    hypertrophic_ECM,
    ossification_ECM,
    pericondrium_ECM
};

// Core (scale factor for lALERT)
float proliferative_CORE = sqrt(1.5);
float hypertrophic_CORE = sqrt(3);
float pericondrium_CORE = 1;
float ossification_CORE = sqrt(3);
float resting_CORE = 1;

float CORE[] = {
    resting_CORE,
    proliferative_CORE,
    hypertrophic_CORE,
    ossification_CORE,
    pericondrium_CORE
};


// Region cell_volume (scale factor for lS)
float proliferative_CV = 1.5;
float hypertrophic_CV = 2.5;
float pericondrium_CV = 1;
float ossification_CV = 2.5;
float resting_CV = 1;

float CV[] = {
    resting_CV,
    proliferative_CV,
    hypertrophic_CV,
    ossification_CV,
    pericondrium_CV
};

// Region line tension (scale factor for eg0)
float proliferative_EG = 0.75;
float hypertrophic_EG = 1;
float pericondrium_EG = 1;
float ossification_EG = 1;
float resting_EG = 1;

float EG[] = {
    resting_EG,
    proliferative_EG,
    hypertrophic_EG,
    ossification_EG,
    pericondrium_EG
};



// Edge threshold
int max_num_edges = 14;
int min_division_n_edges = 12;
int min_exist_n_edges = 6;
//*******************DECLARATIONS*************************
#include "_functions.h"
#include "_policy_fates.h"


void set_solver_threads(int count) {
    num_threads_count = num_threads_division = num_threads_edge_manipulation = count;
    num_threads_force = num_threads_calculation = num_threads_update = count;
}

// Floating-point mechanical and zone parameters accepted by either driver.
struct ModelParameter {
    const char* name;
    float& value;
    float default_value = value;
    float* core = nullptr;
    float default_core = core ? *core : 0;
    void reset() const {
        value = default_value;
        if (core) *core = default_core;
    }
    void set(float next) const {
        value = next;
        if (core) *core = sqrt(next);
    }
};

ModelParameter model_parameters[] = {
    {"k0", k0}, {"kS", kS}, {"lS", lS}, {"kA", kA}, {"e_g0", e_g0},
    {"maxDR", maxDR}, {"lthUP", lthUP}, {"lthDOWN", lthDOWN}, {"lALERT", lALERT},
    {"resting_CV", CV[0], CV[0], &CORE[0]}, {"proliferative_CV", CV[1], CV[1], &CORE[1]},
    {"hypertrophic_CV", CV[2], CV[2], &CORE[2]}, {"ossification_CV", CV[3], CV[3], &CORE[3]},
    {"resting_ECM", ECM[0]}, {"proliferative_ECM", ECM[1]}, {"hypertrophic_ECM", ECM[2]},
    {"ossification_ECM", ECM[3]}, {"resting_T", EG[0]}, {"proliferative_T", EG[1]},
    {"hypertrophic_T", EG[2]}, {"ossification_T", EG[3]}
};

ModelParameter* find_parameter(const std::string& name) {
    for (auto& parameter : model_parameters) if (parameter.name == name) return &parameter;
    return nullptr;
}

class Solver2D {
    bool allocated = false;
    bool decision_mode;
    double clock = 0;
    Tissue2D committed;

    void release_workspace() {
        if (allocated) {
            deallocate();
            force_workspace = {};
        }
        allocated = false;
        release_domain();
    }

    void commit_observation() {
        // Edge manipulation can replace vertices after the CLI sampling point.
        calc_edges();
        calc_areas();
        for (int i = 1; i <= Nc; ++i) if (exist[3][i]) calc_cell_center(i);
        update_shoulder_vertices();
        Tissue2D trial = measure_tissue(clock);
        for (const auto& cell : trial.cells) {
            if (cell.vertices.size() < 3 || !(cell.area > 0) || !std::isfinite(cell.area) ||
                !std::isfinite(cell.target_area) || !std::isfinite(cell.x) || !std::isfinite(cell.y)) {
                throw std::runtime_error("invalid geometry or target area in cell " + std::to_string(cell.id));
            }
            if (training2d.mode) {
                Ring ring;
                for (const auto& p : cell.vertices) ring.push_back({p.x,p.y});
                if (!simple_ring(ring)) throw std::runtime_error("invalid polygon in cell " + std::to_string(cell.id));
            }
            for (const auto& p : cell.vertices) {
                if (!std::isfinite(p.x) || !std::isfinite(p.y) ||
                    !std::isfinite(p.shoulder_x) || !std::isfinite(p.shoulder_y)) {
                    throw std::runtime_error("nonfinite vertex in cell " + std::to_string(cell.id));
                }
            }
        }
        committed = std::move(trial);
    }

public:
    explicit Solver2D(bool decisions) : decision_mode(decisions) {}
    Solver2D(const Solver2D&) = delete;
    Solver2D& operator=(const Solver2D&) = delete;
    ~Solver2D() { release_workspace(); }

    void reset(unsigned seed, const std::string& snapshot = "", int initial_frame = 0) {
        release_workspace();
        committed = {};
        Time = TimeCount = h = 0;
        clock = 0;
        kount = initial_frame;
        placement_seed = seed;
        ddTH = sqrt(pow(lthUP/2,2)+pow(lALERT,2));
        ddTHS = sqrt(pow(lthUP/2,2)+pow(lS+lALERT,2));
        vMAX = 2.*maxDR/h0;
        if (read_and_set_domain() || read_and_set_interface()) {
            throw std::runtime_error("cannot initialize domain and zone interfaces");
        }
        allocate();
        allocated = true;
        if (snapshot.empty()) set_initial_singleCell();
        else if (resume_from_snapshot(snapshot)) throw std::runtime_error("cannot resume snapshot " + snapshot);
        if (training2d.mode) {
            training2d.agents.assign(array_max+1,0);
            training2d.born.assign(array_max+1,0);
            for (int i=1; i<=Nc; ++i) if (exist[3][i]) training2d.birth(i,0);
        }
        if (decision_mode) commit_observation();
    }

    void act(int slot, int action, double angle) {
        const auto& z=training2d.zones.at(training2d.index.at(trained_zone(slot)));
        if (action==0 && z.allow_growth) c_growth[slot]=1;
        else if (action==1 && z.allow_stop) c_growth[slot]=0;
        else if (action==2 && z.allow_division) divide_at_angle(slot,angle,clock);
        else if (action==3 && z.allow_death) { c_dead[slot]=1; c_growth[slot]=0; }
    }

    void training_growth(double dt) {
        const int initial=Nc;
        for (int i=1; i<=initial; ++i) if (exist[3][i] && !c_dead[i]) {
            auto& z=training2d.zones.at(training2d.index.at(trained_zone(i)));
            if (z.allow_growth && c_A0[i] < z.target_size) {
                if (training2d.mode==1 && !c_growth[i]) c_growth[i]=1;
                if (c_growth[i]) c_A0[i]=std::min(z.target_size,c_A0[i]+z.growth_rate*dt*(1-c_A0[i]/z.target_size));
            }
            if (training2d.mode==2) continue;
            double threshold=z.target_size/1.5;
            if (z.allow_death && z.death_rate>0 && c_area[i]>=0.8*threshold &&
                fate_draw(i,"death",-std::expm1(-dt*z.death_rate))) {
                c_dead[i]=1; c_growth[i]=0;
            } else if (z.allow_division && z.division_rate>0 && c_area[i]>threshold &&
                       c[i][0]>=min_division_n_edges && fate_draw(i,"divide",-std::expm1(-dt*z.division_rate))) {
                int other=DIVIDE_CELL(i,0.001,get_subdomain(i));
                c_A0[i]=c_A0[other]=0.5*threshold;
                c_growth[i]=c_growth[other]=0;
                training2d.birth(i,clock);
                training2d.birth(other,clock);
            }
        }
    }

    const Tissue2D& observation() const { return committed; }

    bool terminated() const {
        return std::none_of(committed.cells.begin(), committed.cells.end(),
                            [](const Cell2D& cell) { return !cell.dead; });
    }

    void advance(double deadline, void (*sample)(float) = nullptr) {
        if (!allocated) throw std::runtime_error("reset required before advancing");
        try {
            if (training2d.mode==2) {
                bool any_live=false;
                for (int i=1; i<=Nc; ++i) if (exist[3][i] && !c_dead[i]) any_live=true;
                if (!any_live) { commit_observation(); return; }
            }
            while (clock < deadline && (!decision_mode || !terminated())) {
                if (!training2d.mode) growth_division();
                double start=clock;
                float displacement = integrate_motion(clock, deadline, decision_mode);
                if (sample) sample(displacement);
                if (!std::isfinite(displacement) || displacement > maxDR ||
                    (!decision_mode && displacement == 0)) {
                    throw std::runtime_error("invalid displacement");
                }
                reset_c_shrink();
                vMAX = 2.*maxDR/h;
                for (int i = 1; i <= Nc; ++i) resolve_cross_cell(i);
                edge_manipulation();
                if (training2d.mode) training_growth(clock-start);
                if (decision_mode) commit_observation();
            }
        } catch (...) {
            // The public state is a value snapshot. An unsuccessful trial, including
            // any topology edits, is discarded. Only reset can restart the solver.
            release_workspace();
            clock = committed.time;
            Time = static_cast<float>(clock);
            throw;
        }
    }
};
