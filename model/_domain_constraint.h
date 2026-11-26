/*
void vertex_correction(int vertex_id) {
    float v_x = v[vertex_id][1]; float v_y = v[vertex_id][2];
    //printf("x: %f\n", v_x);
    //printf("y: %f\n", v_y);
    float new_vertex[4] = {0};
    //float polygon_interface_y = 10 * 4;
    int* polygon_boundary;
    if (v_y >= polygon_interface_y) {
       polygon_boundary = e_D_head;
    } else {
       polygon_boundary = e_D_tail;
    }
    float min_dist;
    bool first_assign = true;
    for (int i = 1; i <= polygon_boundary[0]; i++) {
	int j = polygon_boundary[i];
        float u_x = e_normal_unit_vector[j][1];
        float u_y = e_normal_unit_vector[j][2];
        float v_e_D_x = v_D[e_D[j][1]][1];
        float v_e_D_y = v_D[e_D[j][1]][2];
        float dist_j = (v_x - v_e_D_x) * u_x + (v_y - v_e_D_y) * u_y;
        if (first_assign && dist_j < 0) {
	    min_dist = fabs(dist_j);
	    first_assign = false;
	}
        //printf("dist: %f\n", dist_j);
        if (dist_j < 0) {
            if (fabs(dist_j) <= min_dist) {
                //printf("Corrected!");
                new_vertex[0] = 1;
                new_vertex[2] = fabs(dist_j) * u_x + v_x;
                new_vertex[3] = fabs(dist_j) * u_y + v_y;
                min_dist = fabs(dist_j);
            }
        }
    }
    // {1, ...} is corrected, {0, ...} means unchange
    if (new_vertex[0] == 1) {
        v[vertex_id][1] = new_vertex[2];
        v[vertex_id][2] = new_vertex[3];
    }
}
*/
void vertex_correction(int vertex_id) {
    if (training2d.mode && !training2d.zones.empty() && training2d.zones.front().rule) {
        auto wall = wall_distance({v[vertex_id][1],v[vertex_id][2]});
        if (!wall.inside) {
            v[vertex_id][1] = wall.closest[0];
            v[vertex_id][2] = wall.closest[1];
        }
        return;
    }
    float v_x = v[vertex_id][1]; float v_y = v[vertex_id][2];
    //printf("x: %f\n", v_x);
    //printf("y: %f\n", v_y);
    float new_vertex[4] = {0};
    float min_dist;
    bool first_assign = true;
    for (int j = 1; j <= e_D[0][0]; j++) {
        float u_n_x = e_normal_unit_vector[j][1];
        float u_n_y = e_normal_unit_vector[j][2];
	float u_t_x = e_tangent_unit_vector[j][1];
        float u_t_y = e_tangent_unit_vector[j][2];
        float v_e_D_x = v_D[e_D[j][0]][1];
        float v_e_D_y = v_D[e_D[j][0]][2];
	float v_e_D_x_2 = v_D[e_D[j][1]][1];
        float v_e_D_y_2 = v_D[e_D[j][1]][2];
	float len1 = fabs((v_x - v_e_D_x) * u_t_x + (v_y - v_e_D_y) * u_t_y);
        float len2 = fabs((v_x - v_e_D_x_2) * u_t_x + (v_y - v_e_D_y_2) * u_t_y);
	float e_length = torus_distance(v_e_D_x, v_e_D_y, v_e_D_x_2, v_e_D_y_2);
        float dist_j = (v_x - v_e_D_x) * u_n_x + (v_y - v_e_D_y) * u_n_y;
        if (first_assign && dist_j < 0 && len1 <= e_length && len2 <= e_length) {
	    min_dist = fabs(dist_j);
	    first_assign = false;
	}
        //printf("dist: %f\n", dist_j);
        if (dist_j < 0 && len1 <= e_length && len2 <= e_length) {
            if (fabs(dist_j) <= min_dist) {
                //printf("Corrected!");
                new_vertex[0] = 1;
                new_vertex[2] = fabs(dist_j) * u_n_x + v_x;
                new_vertex[3] = fabs(dist_j) * u_n_y + v_y;
                min_dist = fabs(dist_j);
            }
        }
    }
    // {1, ...} is corrected, {0, ...} means unchange
    if (new_vertex[0] == 1) {
        v[vertex_id][1] = new_vertex[2];
        v[vertex_id][2] = new_vertex[3];
    }
}
/*
bool in_domain(float x, float y) {
    //float polygon_interface_y = 10 * 4;
    int* polygon_boundary;
    if (y >= polygon_interface_y) {
       polygon_boundary = e_D_head;
    } else {
       polygon_boundary = e_D_tail;
    }
    for (int i = 1; i <= polygon_boundary[0]; i++) {
	int j = polygon_boundary[i];
        float u_x = e_normal_unit_vector[j][1];
        float u_y = e_normal_unit_vector[j][2];
        float v_e_D_x = v_D[e_D[j][1]][1];
        float v_e_D_y = v_D[e_D[j][1]][2];
        float dist_j = (x - v_e_D_x) * u_x + (y - v_e_D_y) * u_y;
        float min_dist = fabs(dist_j);
        //printf("dist: %f\n", dist_j);
        if (dist_j < 0) {
            return false;
        }
    }
    return true;
}

bool in_domain_some_distance_to_edge(float x, float y) {
    //float polygon_interface_y = 10 * 4;
    int* polygon_boundary;
    if (y >= polygon_interface_y) {
       polygon_boundary = e_D_head;
    } else {
       polygon_boundary = e_D_tail;
    }
    for (int i = 1; i <= polygon_boundary[0]; i++) {
        int j = polygon_boundary[i];
        float u_x = e_normal_unit_vector[j][1];
        float u_y = e_normal_unit_vector[j][2];
        float v_e_D_x = v_D[e_D[j][1]][1];
        float v_e_D_y = v_D[e_D[j][1]][2];
        float dist_j = (x - v_e_D_x) * u_x + (y - v_e_D_y) * u_y;
        //float min_dist = fabs(dist_j);
        //printf("dist: %f\n", dist_j);
	//for points outside region and within pos_dist_thresh to an edge in domain, they should be excluded
	float pos_dist_thresh = 1;
        if (dist_j < pos_dist_thresh) {
            return false;
        }
    }
    return true;
}
*/
bool in_domain_some_distance_to_edge(float x, float y) {
    if (training2d.mode && !training2d.zones.empty() && training2d.zones.front().rule) {
        auto wall = wall_distance({x,y});
        return wall.inside && wall.squared>=1;
    }
    printf("Point %f, %f\n", x, y);
    for (int j = 1; j <= e_D[0][0]; j++) {
        float u_n_x = e_normal_unit_vector[j][1];
        float u_n_y = e_normal_unit_vector[j][2];
        float u_t_x = e_tangent_unit_vector[j][1];
        float u_t_y = e_tangent_unit_vector[j][2];
        float v_e_D_x = v_D[e_D[j][0]][1];
        float v_e_D_y = v_D[e_D[j][0]][2];
        float v_e_D_x_2 = v_D[e_D[j][1]][1];
        float v_e_D_y_2 = v_D[e_D[j][1]][2];
        float len1 = fabs((x - v_e_D_x) * u_t_x + (y - v_e_D_y) * u_t_y);
	float len2 = fabs((x - v_e_D_x_2) * u_t_x +(y - v_e_D_y_2) * u_t_y);
        float dist_j = (x - v_e_D_x) * u_n_x + (y - v_e_D_y) * u_n_y;
        float dist_j_2 = (x - v_e_D_x_2) * u_n_x + (y - v_e_D_y_2) * u_n_y;
	//printf("Point %f, %f\n", x, y);
        //printf("Vector %f, %f\n", (x - v_e_D_x), (y - v_e_D_y));
        //printf("v_n %f, %f\n", u_n_x, u_n_y);
        //printf("dist: %f\n", dist_j);
        //printf("dist: %f\n", dist_j_2);
	//float min_dist = fabs(dist_j);
        //printf("dist: %f\n", dist_j);
	//for points outside region and within pos_dist_thresh to an edge in domain, they should be excluded
	float pos_dist_thresh = 1;
        float e_length = torus_distance(v_e_D_x, v_e_D_y, v_e_D_x_2, v_e_D_y_2);
	//printf("%f\n", e_length);
//	printf("%f\n", len1);
//	printf("%f\n", len2);
	if ((dist_j < pos_dist_thresh || dist_j_2 < pos_dist_thresh) && len1 <= e_length && len2 <= e_length) {
           // printf("Point %f, %f\n", x, y);
	    //printf("Vector %f, %f\n", (x - v_e_D_x), (y - v_e_D_y));
	    //printf("v_n %f, %f\n", u_n_x, u_n_y);
	    //printf("dist: %f\n", dist_j);
	    return false;
        }
    }
    return true;
}
/*
void error_vertex_deletion() {
   for(int i=1; i<=Nv; i++) if(exist[1][i]==1){
       if (!in_domain(v[i][1], v[i][2])) {
            dissolve_vertex(i);
       } 
   }
}
*/
