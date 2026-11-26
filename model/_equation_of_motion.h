#include <omp.h>

void assemble_edge_row(int i, ContactWrites& writes) {
    auto& marks = force_workspace.marked_vertices[i];
    for (int j : marks) edge_vertex_int[i][j] = 0;
    marks.clear();
    int cell = e_cells[i][1];
    int a = e[i][1], b = e[i][2];
    if (e_sign_cell(i, cell) < 0) std::swap(a, b);
    auto& due = force_workspace.edge_due[i];
    for (size_t block = 0; block < due.size(); ++block) {
        if (contact_now() < due[block]) continue;
        float earliest = std::numeric_limits<float>::infinity();
        int first = block * pair_block_size + 1, last = std::min(first + pair_block_size - 1, Nv);
        for (int j = first; j <= last; ++j) if (exist[1][j] && v_cells[j][1] != cell) {
            if (contact_now() >= ve_next_time[i][j]) assemble_vertex_edge(i, j, a, b, cell, writes);
            earliest = std::min(earliest, ve_next_time[i][j]);
        }
        due[block] = earliest;
    }
}

void assemble_shoulder_row(int i, ContactWrites& writes) {
    int cell = v_cells[i][1];
    auto& due = force_workspace.shoulder_due[i];
    for (size_t block = i / pair_block_size; block < due.size(); ++block) {
        if (contact_now() < due[block]) continue;
        float earliest = std::numeric_limits<float>::infinity();
        int first = block * pair_block_size + 1, last = std::min(first + pair_block_size - 1, Nv);
        for (int j = std::max(i + 1, first); j <= last; ++j) if (exist[1][j] && v_cells[j][1] != cell) {
            if (contact_now() >= vv_next_time[i][j]) assemble_shoulder(i, j, cell, v_cells[j][1], writes);
            earliest = std::min(earliest, vv_next_time[i][j]);
        }
        due[block] = earliest;
    }
}

void apply_core_contacts(const ContactWrites& writes) {
    if (training2d.mode == 2) return;
    for (const auto& pair : writes.core_contacts) {
        int a=pair.first, b=pair.second;
        c_growth[a] = c_growth[b] = 0;
        if (!c_shrink[a] && c_A0[a] >= 0.5*CV[force_workspace.zone[a]]) {
            c_A0[a] *= 0.85;
            c_shrink[a] = 1;
        }
        if (!c_shrink[b] && c_A0[b] >= 0.5*CV[force_workspace.zone[b]]) {
            c_A0[b] /= 0.85;
            c_shrink[b] = 1;
        }
    }
}

float calc_forces() {
    auto& work = force_workspace;
    work.prepare();
    float maximum = 0;
    #pragma omp parallel num_threads(num_threads_force) reduction(max:maximum)
    {
        int thread = omp_get_thread_num();
        int team = omp_get_num_threads();
        auto& writes = work.contact[thread];
        writes.prepare(Nv);
        #pragma omp for schedule(static)
        for (int i = 1; i <= Nc; ++i) if (exist[3][i]) {
            work.zone[i] = get_subdomain(i);
            assemble_area_force(i, work);
        }
        #pragma omp for schedule(static) nowait
        for (int i = 1; i <= Ne; ++i) if (exist[2][i]) assemble_length_force(i, work);

        // Static contiguous row ranges make thread-index order match the
        // reference's outer loop. Inside each row, j always increases.
        #pragma omp for schedule(static)
        for (int i = 1; i <= Ne; ++i) if (exist[2][i]) assemble_edge_row(i, writes);
        #pragma omp for schedule(static)
        for (int i = 1; i <= Nv; ++i) if (exist[1][i]) assemble_shoulder_row(i, writes);

        #pragma omp single nowait
        {
            wA = wL = wS = 0;
            for (int i = 1; i <= Nc; ++i) if (exist[3][i]) wA += work.area_energy[i];
            for (int i = 1; i <= Ne; ++i) if (exist[2][i]) wL += work.length_energy[i];
            for (int t = 0; t < team; ++t) apply_core_contacts(work.contact[t]);
        }
        #pragma omp for schedule(static)
        for (int i = 1; i <= Nv; ++i) if (exist[1][i]) {
            float* force = v_F[i];
            for (const auto& term : work.area_terms[i]) term.add(force);
            for (int edge : work.incident_edges[i]) work.length_terms[edge][e[edge][2] == i].add(force);
            for (int t = 0; t < team; ++t)
                for (const auto& delta : work.contact[t].edge_vertex[i]) delta.add(force);
            for (int t = 0; t < team; ++t)
                for (const auto& delta : work.contact[t].shoulder[i]) delta.add(force);
            float magnitude = sqrt(force[1]*force[1] + force[2]*force[2]);
            if (magnitude > maximum) maximum = magnitude;
        }
    }
    return maximum;
}
void calc_edges(){
#pragma omp parallel for num_threads(num_threads_calculation)
    for(int i=1; i<=Ne; i++) if(exist[2][i]==1){
        e_length[i]=edge_length(i);
        edge_center(i);
        tangent_and_normal(i);
    }
}
//****************************************************************************
void calc_areas(){
#pragma omp parallel for num_threads(num_threads_calculation)
    for(int i=1; i<=Nc; i++) if(exist[3][i]==1) {
        c_area[i]=CellArea_new(i);
	//if (c_area[i] < 0.1) {
	//    dissolve_cell(i);
	//}
    } 
}
//****************************************************************************
void update_edges(){
#pragma omp parallel for num_threads(num_threads_update)
    for(int i=1; i<=Ne; i++) if(exist[2][i]==1) {
        float elen=edge_length(i);
        e_dl[i]=(elen-e_length[i]);
        e_length[i]=elen;
        edge_center(i);
        tangent_and_normal(i);
    }
}
//****************************************************************************
void update_areas(){
#pragma omp parallel for num_threads(num_threads_update)
    for(int i=1; i<=Nc; i++) if(exist[3][i]==1) c_area[i]=CellArea_new(i);
}
//****************************************************************************
// The centroid of ECM will not be the same as the cell, and the shape will change, but this simplifies the calculation of force
void update_shoulder_vertices(){
#pragma omp parallel for num_threads(num_threads_update)
    for(int i=1; i<=Nv; i++) if(exist[1][i]==1){
        int e1=v_edges[i][1];
        int e2=v_edges[i][2];
	float cell_id = v_cells[i][1];
	float ss = ECM[get_subdomain(cell_id)];
        float nrm=sqrt(pow(e_n[e1][1]+e_n[e2][1],2)+pow(e_n[e1][2]+e_n[e2][2],2));
        vS[i][1]=v[i][1]+lS*ss*(e_n[e1][1]+e_n[e2][1])/nrm;
        vS[i][2]=v[i][2]+lS*ss*(e_n[e1][2]+e_n[e2][2])/nrm;
        /*
	if (v[i][1] == 0 && v[i][2] == 0) {
	    int e1 = v_edges[i][1];
	    int e2 = v_edges[i][2];
	    int v11 = e[e1][1];
	    int v12 = e[e1][2];
	    int v21 = e[e2][1];
	    int v22 = e[e2][2];
	    int vend1;
	    int vend2;
	    if (v11 != i) {
	        vend1 = v11;
	    } else {
	    	vend1 = v12;
	    }
	    if (v21 != i) {
	    	vend2 = v21;
	    } else {
	    	vend2 = v22;
	    }
	    make_edge(vend1, vend2, 0);
	    dissolve_edge(e1, 1);
	    dissolve_edge(e2, 1);
	    dissolve_vertex(i);
	}
       */	
        // You can preserve the centroid using these statements
        //vS[i][1]=v[i][1]+lS*(e_n[e1][1]+e_n[e2][1]);
        //vS[i][2]=v[i][2]+lS*(e_n[e1][2]+e_n[e2][2]);
    }
}
//****************************************************************************
// The CLI retains its original single-precision step and clock. Decision-mode
// deadlines use double precision, including the final displacement multiplier.
float integrate_motion(double& clock, double deadline, bool decision_mode) {
    calc_edges();
    calc_areas();
    float FMAX = calc_forces();
    h = (maxDR * 0.9) / FMAX;
    double dt = h;
    if (decision_mode) {
        if (!std::isfinite(FMAX)) throw std::runtime_error("nonfinite force");
        if (FMAX == 0) dt = h0;
        dt = std::min(dt, deadline - clock);
        if (!(dt > 0) || !std::isfinite(dt) || clock + dt == clock) {
            throw std::runtime_error("invalid adaptive step");
        }
        h = static_cast<float>(dt);
    } else {
        printf("Time step chosen is %f", h);
    }

    float displacement = FMAX * h;
    contact_travel += static_cast<double>(FMAX) * dt;
    #pragma omp parallel for num_threads(num_threads_update)
    for (int i = 1; i <= Nv; ++i) if (exist[1][i]) {
        if (!c_dead[v_cells[i][1]]) {
            for (int j = 1; j <= 2; ++j) {
                if (decision_mode) v[i][j] += dt * v_F[i][j];
                else v[i][j] += h * v_F[i][j];
            }
        }
        vertex_correction(i);
        v_F[i][1] = v_F[i][2] = 0;
    }
    // A center can read a vertex that is still awaiting its second correction.
    // Keep this pass in dev order for walls whose second projection moves it.
    for (int i = 1; i <= Nv; ++i) if (exist[1][i]) {
        vertex_correction(i);
        // Preserve the legacy traversal and its cached centers for CLI frames.
        if (c_vertices[i][0] != 0) calc_cell_center(i);
    }
    update_edges();
    update_areas();
    update_shoulder_vertices();

    if (decision_mode) {
        clock = dt == deadline - clock ? deadline : clock + dt;
        Time = static_cast<float>(clock);
    } else {
        Time += h;
        clock = Time;
    }
    return displacement;
}
