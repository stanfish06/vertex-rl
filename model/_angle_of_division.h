#include <cmath>

float compute_division_angle(int i) {
    float c_x = c_center[i][1];

    float b_l = boundary_x[1];
    float b_r = boundary_x[2];
	
    float n_l_x = e_normal_unit_vector[lateral_normal[0]][1];
    float n_l_y = e_normal_unit_vector[lateral_normal[0]][2];
    float n_r_x = -e_normal_unit_vector[lateral_normal[1]][1];
    float n_r_y = -e_normal_unit_vector[lateral_normal[1]][2];
    //printf("n_l_x: %f", n_l_x);
    //printf("n_r_x: %f", n_r_x);
    float n_d_x = (c_x - b_l) / (b_r - b_l) * n_l_x + (b_r - c_x) / (b_r - b_l) * n_r_x;
    float n_d_y = (c_x - b_l) / (b_r - b_l) * n_l_y + (b_r - c_x) / (b_r - b_l) * n_r_y;
    // Angle of vector orthogonal to normal vector
    float a = atan(n_d_y / n_d_x);
    if (a < 0) a = a + M_PI;
    return a;
}
