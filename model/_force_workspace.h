// Intermediate force terms retain the precision of the legacy expressions.
// Each addition rounds back to the float vertex accumulator in reference order.
struct Force2 {
    double x, y;
    void add(float* target) const {
        target[1] += x;
        target[2] += y;
    }
};

struct ContactWrites {
    std::vector<std::vector<Force2>> edge_vertex, shoulder;
    std::vector<std::pair<int, int>> core_contacts;

    void prepare(int vertices) {
        edge_vertex.resize(vertices + 1);
        shoulder.resize(vertices + 1);
        for (auto& row : edge_vertex) row.clear();
        for (auto& row : shoulder) row.clear();
        core_contacts.clear();
    }
};

// A block may be skipped only until its earliest dev next-time entry. Topology
// edits invalidate the minima because cell ownership and reused IDs may change.
constexpr int pair_block_size = 32;

struct ForceWorkspace {
    std::vector<std::array<Force2, 2>> area_terms, length_terms;
    std::vector<float> area_energy, length_energy;
    std::vector<int> zone;
    std::vector<std::vector<int>> incident_edges, marked_vertices;
    std::vector<std::vector<float>> edge_due, shoulder_due;
    std::vector<ContactWrites> contact;

    void prepare() {
        contact.resize(num_threads_force);
        if (!force_topology_changed) return;
        force_topology_changed = false;
        area_terms.assign(Nv + 1, {});
        length_terms.resize(Ne + 1);
        area_energy.resize(Nc + 1);
        length_energy.resize(Ne + 1);
        zone.resize(Nc + 1);
        incident_edges.resize(Nv + 1);
        for (auto& row : incident_edges) row.clear();
        for (int i = 1; i <= Ne; ++i) if (exist[2][i]) {
            incident_edges[e[i][1]].push_back(i);
            incident_edges[e[i][2]].push_back(i);
        }
        marked_vertices.resize(Ne + 1);
        size_t blocks = (Nv + pair_block_size - 1) / pair_block_size;
        edge_due.resize(Ne + 1);
        shoulder_due.resize(Nv + 1);
        for (auto& row : edge_due) row.assign(blocks, -std::numeric_limits<float>::infinity());
        for (auto& row : shoulder_due) row.assign(blocks, -std::numeric_limits<float>::infinity());
    }
};

ForceWorkspace force_workspace;
