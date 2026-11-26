//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
float vv_steric_shoulder_force(int i, int j){
    int ei1=v_edges[i][1];
    int ei2=v_edges[i][2];
    int ej1=v_edges[j][1];
    int ej2=v_edges[j][2];
    float *dxdydz = new float[3];
    int cell1=v_cells[i][1];
    int cell2=v_cells[j][1];
    float cell1_A = c_area[cell1];
    float cell2_A = c_area[cell2];
    float cell1_CV = CV[get_subdomain(cell1)];
    float cell2_CV = CV[get_subdomain(cell2)];
    float ss_c1 = ECM[get_subdomain(cell1)];
    float ss_c2 = ECM[get_subdomain(cell2)]; 
    if( edge_vertex_int[j][ei1]==0 && edge_vertex_int[j][ei2]==0){
        float dd=vv_distance(i,j);
        vv_next_time[i][j]=Time+(dd-lS)/vMAX;
        vv_next_time[j][i]=vv_next_time[i][j];
        if( dd<=lS*ss_c1 || dd<=lS*ss_c2 ){
            //i (reference vertex)
            float vix=v[i][1];
            float viy=v[i][2];
            //j
            //dxdydz[1]=0; dxdydz[2]=0;
            //torus_dx_dy_dz(dxdydz,j,i);
            float vjx=v[j][1];
            float vjy=v[j][2];
	    // either one edges' normal would be fine
	    //printf("centroid is x: %f | y: %f\n", cell2_c_x, cell2_c_y);
	    //printf("normal is x: %f | y: %f\n", e_n[ej1][1], e_n[ej1][2]);
	    //printf("edge vertex is x: %f | y: %f", vjx, vjy);
	    if (dd<=lS*ss_c2 && cell1_A >= 0.3 * cell1_CV) {
                // If dd > IS, there will be weak attraction
                //i
                v_F[i][1]+=(2*kS*(vix - vjx)*(lS*ss_c2 - dd))/dd;
                v_F[i][2]+=(2*kS*(lS*ss_c2 - dd)*(viy - vjy))/dd;
		//v_F[i][1]+=(2*kS*(vix - vjx)*(lS*ss_c2 - dd))/dd;
                //v_F[i][2]+=(2*kS*(lS*ss_c2 - dd)*(viy - vjy))/dd;
	    }
            //j
	    if (dd<=lS*ss_c1 && cell2_A >= 0.3 * cell2_CV) {
                v_F[j][1]+=(2*kS*(vix - vjx)*(-lS*ss_c1 + dd))/dd;
                v_F[j][2]+=(2*kS*(-lS*ss_c1 + dd)*(viy - vjy))/dd;
		//v_F[j][1]+=(2*kS*(vix - vjx))/dd;
                //v_F[j][2]+=(2*kS*(viy - vjy))/dd;
	    }
        }
    }
    
    delete []dxdydz;
    
    return 0;
}
//****************************************************************************
