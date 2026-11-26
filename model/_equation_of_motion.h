//****************************************************************************
//****************************************************************************
//***********************EQUATION OF MOTION***********************************
//****************************************************************************
//****************************************************************************
#include <omp.h>
float calc_forces(){    
float Fmax=0; 
#pragma omp parallel num_threads(num_threads_force)
    {
	int thread_id = omp_get_thread_num();
    	wA=0;
    	wL=0;
    	wS=0;
//	if (thread_id == 0) {
            //AREA COMPRESSIBILITY
   	    #pragma omp for schedule(dynamic, 1)
       	    for(int i=1; i<=Nc; i++) if(exist[3][i]==1) wA+=c_AreaCompressibility_force_New(i);
//	} else if (thread_id == 1) {
    	    //LINE ENERGY
   	    #pragma omp for schedule(dynamic, 1)
    	    for(int i=1; i<=Ne; i++) if(exist[2][i]==1) wL += e_length_force(i);
//    	} else if (thread_id == 2) {
    	    //VERTEX-EDGE STERIC INTERACTION & SHOULDER INTERACTION
    	    #pragma omp for schedule(dynamic, 1)
    	    for(int i=1; i<=Ne; i++) if(exist[2][i]==1){
                int eCELL=e_cells[i][1];
                for(int j=1; j<=Nv; j++) if(exist[1][j]==1){
            	    edge_vertex_int[j][i]=0;
            	    int vCELL=v_cells[j][1];
            	    if(vCELL!=eCELL && Time>=ve_next_time[i][j]) wS+=ve_steric_force(i,j);
       	        }
    	    }
//   	}
//	} else if (thread_id == 3) {

	#pragma omp barrier    
	#pragma omp parallel for num_threads(num_threads_force)
    	for(int i=1; i<=Nv; i++) if(exist[1][i]==1){
            int vCELL1=v_cells[i][1];
            for(int j=i+1; j<=Nv; j++) if(exist[1][j]==1){
                int vCELL2=v_cells[j][1];
                if(vCELL1!=vCELL2 && Time>=vv_next_time[i][j]) wS+=vv_steric_shoulder_force(i,j);
            }
        }
//	} 
	//else if (thread_id == 4) {
    	    //MAX FORCE
	   // #pragma omp barrier
    	   // #pragma omp parallel for num_threads(num_threads_force)
    	  //  for(int i = 1; i <= Nv; i++) if(exist[1][i]==1){
            //    float foRce=sqrt(v_F[i][1]*v_F[i][1]+v_F[i][2]*v_F[i][2]);
              //  if(foRce>Fmax) Fmax=foRce;
            //}
	//}
    
	//MAX FORCE
    	#pragma omp barrier
    	#pragma omp parallel for num_threads(num_threads_force)
    	for(int i = 1; i <= Nv; i++) if(exist[1][i]==1){
            float foRce=sqrt(v_F[i][1]*v_F[i][1]+v_F[i][2]*v_F[i][2]);
            if(foRce>Fmax) Fmax=foRce;
        }
    }

    return Fmax;
}
//****************************************************************************
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
    #pragma omp parallel for num_threads(num_threads_update)
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
