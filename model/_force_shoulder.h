void assemble_shoulder(int i, int j, int cell1, int cell2, ContactWrites& writes) {
    if (edge_vertex_int[v_edges[i][1]][j] || edge_vertex_int[v_edges[i][2]][j]) return;
    float dd = vv_distance(i,j);
    int zone1 = force_workspace.zone[cell1], zone2 = force_workspace.zone[cell2];
    vv_next_time[i][j] = contact_due(dd, lS*std::max(ECM[zone1], ECM[zone2]), lS);
    float ss_c1 = ECM[zone1], ss_c2 = ECM[zone2];
    float vix=v[i][1], viy=v[i][2], vjx=v[j][1], vjy=v[j][2];
    if (dd <= lS*ss_c2 && c_area[cell1] >= 0.3*CV[zone1]) {
        writes.shoulder[i].push_back({(2*kS*(vix-vjx)*(lS*ss_c2-dd))/dd,
                                     (2*kS*(lS*ss_c2-dd)*(viy-vjy))/dd});
    }
    if (dd <= lS*ss_c1 && c_area[cell2] >= 0.3*CV[zone2]) {
        writes.shoulder[j].push_back({(2*kS*(vix-vjx)*(-lS*ss_c1+dd))/dd,
                                     (2*kS*(-lS*ss_c1+dd)*(viy-vjy))/dd});
    }
}
