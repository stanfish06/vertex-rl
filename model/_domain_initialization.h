void compute_normal_vector() {
    for (int i = 1; i <= e_D[0][0]; i++) {
        float v_e_D_1_x = v_D[e_D[i][0]][1];
        float v_e_D_1_y = v_D[e_D[i][0]][2];
        float v_e_D_2_x = v_D[e_D[i][1]][1];
        float v_e_D_2_y = v_D[e_D[i][1]][2];

        float length_e = sqrt(pow(v_e_D_2_x - v_e_D_1_x, 2) + pow(v_e_D_2_y - v_e_D_1_y, 2));
        float u_e_D_x = (v_e_D_2_y - v_e_D_1_y) / length_e;
        float u_e_D_y = (v_e_D_1_x - v_e_D_2_x) / length_e;

        float v_e_D_mid_x = (v_e_D_1_x + v_e_D_2_x) / 2;
        float v_e_D_mid_y = (v_e_D_1_y + v_e_D_2_y) / 2;
        
        float max_dist_to_another_v_D = 0;
        int most_distant_v_D = 0;
        float dist_j = 0;
        for (int j = 1; j <= v_D[0][0]; j++) {
            if (j != e_D[i][0] && j != e_D[i][1]) {
                dist_j = sqrt(pow(v_D[j][1] - v_e_D_mid_x, 2) + pow(v_D[j][2] - v_e_D_mid_y, 2));
                if (dist_j > max_dist_to_another_v_D) {
                    max_dist_to_another_v_D = dist_j;
                    most_distant_v_D = j;
                }
            }
        }

        float perturb_dist = 1.5;
        v_e_D_mid_x = v_e_D_mid_x + u_e_D_x * perturb_dist;
        v_e_D_mid_y = v_e_D_mid_y + u_e_D_y * perturb_dist;
        dist_j = sqrt(pow(v_D[most_distant_v_D][1] - v_e_D_mid_x, 2) + pow(v_D[most_distant_v_D][2] - v_e_D_mid_y, 2));
        if (dist_j > max_dist_to_another_v_D) {
            u_e_D_x = -u_e_D_x;
            u_e_D_y = -u_e_D_y;
        }
        e_normal_unit_vector[i][1] = u_e_D_x;
        e_normal_unit_vector[i][2] = u_e_D_y;
    }
}

void print_normal_vector() {
    for (int i = 1; i <= e_normal_unit_vector[0][0]; i++) {
        printf("x: %f y: %f\n", e_normal_unit_vector[i][1], e_normal_unit_vector[i][2]);
    }
}