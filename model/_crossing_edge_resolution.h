int orientation(int v1, int v2, int v3) {
    float x_v1 = v[v1][1];
    float y_v1 = v[v1][2];
    float x_v2 = v[v2][1];
    float y_v2 = v[v2][2];
    float x_v3 = v[v3][1];
    float y_v3 = v[v3][2];

    float cross_prod = (x_v3 - x_v1) * (y_v2 - y_v1) - (x_v2 - x_v1) * (y_v3 - y_v1);
    if (cross_prod > 0) return 1;
    else return -1;
}

bool is_cross(int e1, int e2) {
    int v11 = e[e1][1];
    int v12 = e[e1][2];
    int v21 = e[e2][1];
    int v22 = e[e2][2];
    if ((v11 == v21) || (v11 == v22) || (v12 == v21) || (v12 == v22)) {
        return false;
    }
    
    int orient_e1_v21 = orientation(v11, v12, v21);
    int orient_e1_v22 = orientation(v11, v12, v22);
    int orient_e2_v11 = orientation(v21, v22, v11);
    int orient_e2_v12 = orientation(v21, v22, v12);
    /*
    printf("e1:\n");
        printf("x11: %f x12: %f\n", v[v11][1], v[v12][1]);
        printf("y11: %f y12: %f\n", v[v11][2], v[v12][2]);
        printf("e2:\n");
        printf("x21: %f x22: %f\n", v[v21][1], v[v22][1]);
        printf("y21: %f y22: %f\n", v[v21][2], v[v22][2]);
        printf("%d, %d\n", orient_e1_v21, orient_e1_v22);
        printf("%d, %d\n", orient_e2_v11, orient_e2_v12);
    */
    if ((orient_e1_v21 != orient_e1_v22) && (orient_e2_v11 != orient_e2_v12)) {
        /*
	printf("e1:\n");
        printf("x11: %f x12: %f\n", v[v11][1], v[v12][1]);
        printf("y11: %f y12: %f\n", v[v11][2], v[v12][2]);
	printf("e2:\n");
	printf("x21: %f x22: %f\n", v[v21][1], v[v22][1]);
        printf("y21: %f y22: %f\n", v[v21][2], v[v22][2]);
	printf("%d, %d\n", orient_e1_v21, orient_e1_v22);
	printf("%d, %d\n", orient_e2_v11, orient_e2_v12);
	*/
        return true;
    } else {
        return false;
    }
} 

// If a cell gets a cross, the flipped part will be eliminated
void resolve_cross_cell(int cell_id) {
    printf("Start resolving cell %d\n", cell_id);
    /* 
    int num_crosses = 0;
    if (c[cell_id][0] != 0) {
    	for (int i = 1; i <= c[cell_id][0]; i++) {
            for (int j = 1; j <= c[cell_id][0]; j++) {
                if (is_cross(abs(c[cell_id][i]),abs(c[cell_id][j]))) {
                    num_crosses++;
                } 
            }
        }
        printf("Number of crosses: %d", num_crosses);
        num_crosses = num_crosses / 2;
        if (num_crosses > 0) {
           printf("Dissolve the cell");
           dissolve_cell(cell_id);
           delete_cell(cell_id);
           return;
        }
    }
    return;
    */

    bool all_resolved = false;
//rep:    
    if (c[cell_id][0] > 3) {
	//all_resolved = true;
    	for (int i = 1; i <= c[cell_id][0]; i++) {
            for (int j = 1; j <= c[cell_id][0]; j++) {
	        if (is_cross(abs(c[cell_id][i]),abs(c[cell_id][j]))) {
	            printf("Find crossing\n");
		    //printf("Between %d\n", i);
		    //printf("and %d\n", j);
		    e_length[abs(c[cell_id][i])] = lthDOWN / 2;
		    e_length[abs(c[cell_id][j])] = lthDOWN / 2;
		    /*
	            int vi2 = e[i][2];
		    int e1 = v_edges[vi2][1];
		    int e2 = v_edges[vi2][2];
		    int e_next = 0;
		    if (e1 != i) e_next = e1;
		    else e_next = e2;
		    float e_t_x = e_t[i][1]; 
	            float e_t_y = e_t[i][2];
		    float en_t_x = e_t[e_next][1];
		    float en_t_y = e_t[e_next][2];
		    float cross_prod = e_t_x * en_t_y - en_t_x * e_t_y;
		    int pre_e = i;
		    int last_e = j;
		    int current_e = e_next;
		    if (cross_prod < 0) {
		   	pre_e = j;
		   	last_e = i;
		   	int vj2 = e[j][2];
		   	e1 = v_edges[vj2][1];
                        e2 = v_edges[vj2][2];
                        if (e1 != j) current_e = e1;
                        else current_e = e2;
		    }
		    int last_merged_e = 0;
		    while (last_merged_e != last_e) {
			printf("Merge edge %d\n", pre_e);
		        int vce2 = e[current_e][2];
                        e1 = v_edges[vce2][1];
                        e2 = v_edges[vce2][2];
			//EDGE_MERGE(pre_e);
			e_length[pre_e] = lthDOWN / 2;
			last_merged_e = pre_e;
                        pre_e = current_e;
			if (e1 != current_e) current_e = e1;
                        else current_e = e2;
		    }  
	            //all_resolved = false;
		    //goto rep;
		    */ 
		    return; 
	        }
	    }
        }
    }
}
