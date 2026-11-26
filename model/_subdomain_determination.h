/*
int get_subdomain(int cell_id) {
   float c_x = c_center[cell_id][1];
   float c_y = c_center[cell_id][2];
   int* polygon_boundary;
   if (c_y >= polygon_interface_y) {
      polygon_boundary = e_D_head;
   } else {
      polygon_boundary = e_D_tail;
   }
   float min_dist;
   for (int i = 1; i <= polygon_boundary[0]; i++) {
       int j = polygon_boundary[i];
       float u_x = e_normal_unit_vector[j][1];
       float u_y = e_normal_unit_vector[j][2];
       float v_e_D_x = v_D[e_D[j][1]][1];
       float v_e_D_y = v_D[e_D[j][1]][2];
       float dist_j = (c_x - v_e_D_x) * u_x + (c_y - v_e_D_y) * u_y;
       if (i == 1) {
           min_dist = dist_j;
       }
       if (dist_j < min_dist) {
           min_dist = dist_j;
       }
   }
   if (min_dist <= pericondrium_ECM) {
       return 4;
   }
   if (c_y >= resting_zone_interface[0] && c_y <= resting_zone_interface[1]) {
       //printf("resting cells!");
       return 0;
   } else if (c_y >= proliferative_zone_interface[0] && c_y <= proliferative_zone_interface[1]) {
       return 1;
   } else if (c_y >= hypertrophic_zone_interface[0] && c_y <= hypertrophic_zone_interface[1]) {
       return 2;
   } else if (c_y >= ossification_interface[0] && c_y <= ossification_interface[1]) {
       return 3;
   } else {
       return 1; // This shouldn't happen, if it does, check previous conditions
   }
   return 1; // shouldn't be used
}
*/
int zone_at(float c_x, float c_y) {
   bool first_assign = true;
   float min_dist;
   int cell_subdomain = -1;
   for (int i = 1; i <= v_interface[0][0]; i++) {
       float dist = (c_x - v_interface[i][1]) * e_normal_unit_vector_interface[i][1] + (c_y - v_interface[i][2]) * e_normal_unit_vector_interface[i][2];
       if (dist > 0 && first_assign) {
           first_assign = false;
	   min_dist = dist;
	   cell_subdomain = i - 1;
       } else if (dist > 0 && !first_assign && dist < min_dist) {
      	   min_dist = dist;
	   cell_subdomain = i - 1; 
       }
   }
   if (first_assign) {
       cell_subdomain = 3;
   }
   return cell_subdomain;
}

int get_subdomain(int cell_id) {
    if (training2d.mode && !training2d.zones.empty() && training2d.zones.front().rule) {
        PlanarPoint center{c_center[cell_id][1],c_center[cell_id][2]};
        for (const auto& z : training2d.zones) if (zone_matches(z,center)) {
            // Custom zone IDs outside the legacy material table use its base material.
            return z.id < 5 ? z.id : 0;
        }
        // The serial observation check reports an uncovered center.
        return 0;
    }
    return zone_at(c_center[cell_id][1], c_center[cell_id][2]);
}
