void assemble_vertex_edge(int i, int j, int v1, int v2, int cID, ContactWrites& writes) {
    int cell_j = v_cells[j][1];
    int zone_i = force_workspace.zone[cID], zone_j = force_workspace.zone[cell_j];
    float celli_A = c_area[cID], cellj_A = c_area[cell_j];
    float celli_CV = CV[zone_i], cellj_CV = CV[zone_j];
    float ss_lA_ci = CORE[zone_i], ss_lA_cj = CORE[zone_j];
    float ss_ci = ECM[zone_i], ss_cj = ECM[zone_j];

    float dd = torus_distance(v[j][1], v[j][2], e_center[i][1], e_center[i][2]);
    ve_next_time[i][j] = Time+(dd-ddTH)/vMAX;
    float cutoff = sqrt(pow(lthUP/2,2)+pow(lS*std::max(ss_ci,ss_cj)+lALERT*std::max(ss_lA_ci,ss_lA_cj),2));
    if (!(dd < cutoff)) return;
    float v1x=v[v1][1], v1y=v[v1][2], v2x=v[v2][1], v2y=v[v2][2];
    float vx=v[j][1], vy=v[j][2];
    float Rx=vx-v1x, Ry=vy-v1y;
    float RdotN=e_n[i][1]*Rx+e_n[i][2]*Ry;
    float RdotT=e_t[i][1]*Rx+e_t[i][2]*Ry;
    if (!(RdotT > 0 && RdotT < e_length[i])) return;

    // These products retain the dev association and precision. In particular,
    // the end projection contains double-precision pow terms.
    float dx = v1x-v2x, dy = v1y-v2y;
    float start_projection = dx*(v2x-vx) + dy*(v2y-vy);
    auto end_projection = [&] {
        return pow(v1x,2) + pow(v1y,2) + v2x*vx - v1x*(v2x+vx) + v2y*vy - v1y*(v2y+vy);
    };
    if (RdotN < 0 && RdotN > -std::max(ss_lA_ci,ss_lA_cj)*lALERT) {
        float l2=pow(e_length[i],2), l4=l2*l2;
        float cross = -(v2y*vx) + v1y*(-v2x+vx) + v1x*(v2y-vy) + v2x*vy;
        float opposite = v1y*(v2x-vx) + v2y*vx - v2x*vy + v1x*(-v2y+vy);
        double end = end_projection();
        writes.edge_vertex[v1].push_back({-(k0*(2*dy*start_projection*cross)/l4),
                                         -(k0*(-2*dx*start_projection*cross)/l4)});
        writes.edge_vertex[v2].push_back({-(k0*(2*dy*opposite*end)/l4),
                                         -(k0*(2*dx*cross*end)/l4)});
        writes.edge_vertex[j].push_back({-(k0*(2*(-v1y+v2y)*opposite)/l2),
                                        -(k0*(2*dx*opposite)/l2)});
        writes.core_contacts.emplace_back(cID, cell_j);
    }
    if (RdotN > 0 && (RdotN < lS*ss_ci || RdotN < lS*ss_cj)) {
        float li=e_length[i], l4=pow(li,4);
        auto shell_gap = [&](float scale) {
            return lS*scale*li - v2y*vx + v1y*(-v2x+vx) + v1x*(v2y-vy) + v2x*vy;
        };
        if (RdotN < lS*ss_cj && celli_A >= 0.3*celli_CV) {
            float gap = shell_gap(ss_cj);
            double end = end_projection();
            writes.edge_vertex[v1].push_back({-(2*kS*dy*start_projection*gap/l4),
                                             -(-2*kS*dx*start_projection*gap/l4)});
            writes.edge_vertex[v2].push_back({-(-2*kS*dy*gap*end/l4),
                                             -(2*kS*dx*gap*end/l4)});
        }
        if (RdotN < lS*ss_ci && cellj_A >= 0.3*cellj_CV) {
            float gap = shell_gap(ss_ci);
            double squared = pow(dx,2) + pow(dy,2);
            writes.edge_vertex[j].push_back({-(2*kS*dy*gap/squared),
                                            -(2*kS*(-v1x+v2x)*gap/squared)});
        }
        edge_vertex_int[i][j] = 1;
        force_workspace.marked_vertices[i].push_back(j);
    }
}
