//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
float ve_steric_force(int i, int j){
    
    //determining v1 & v2
    int v1=e[i][1], v2=e[i][2];
    int cID=e_cells[i][1];
    int sgn=e_sign_cell(i,cID);
    //printf("edge sign is %d", sgn);
    if(sgn<0){
        //printf("edge sign is %d", sgn);
        int vv=v1;
        v1=v2;
        v2=vv;
    }

    int cell_j = v_cells[j][1];

    float celli_A = c_area[cID];
    float cellj_A = c_area[cell_j];
    float celli_CV = CV[get_subdomain(cID)];
    float cellj_CV = CV[get_subdomain(cell_j)];

    float ss_lA_ci = CORE[get_subdomain(cID)];
    float ss_lA_cj = CORE[get_subdomain(cell_j)];

    float ss_ci = ECM[get_subdomain(cID)];
    float ss_cj = ECM[get_subdomain(cell_j)];
    
    float *dxdydz = new float[3];
    
    float dd=torus_distance(v[j][1],v[j][2],e_center[i][1],e_center[i][2]);
    ve_next_time[i][j]=Time+(dd-ddTH)/vMAX;
    ddTHS=sqrt(pow(lthUP/2,2)+pow(lS*std::max(ss_ci, ss_cj)+lALERT*std::max(ss_lA_ci, ss_lA_cj),2));
    if(dd<ddTHS) {
            
        //v1 (reference vertex)
        float v1x=v[v1][1];
        float v1y=v[v1][2];
        
        //v2
        dxdydz[1]=0; dxdydz[2]=0;
        //torus_dx_dy_dz(dxdydz,v2,v1);
        float v2x=v[v2][1]+dxdydz[1];
        float v2y=v[v2][2]+dxdydz[2];
        
        //j
        dxdydz[1]=0; dxdydz[2]=0;
        //torus_dx_dy_dz(dxdydz,j,v1);
        float vx=v[j][1]+dxdydz[1];
        float vy=v[j][2]+dxdydz[2];
        
        //vector R
        float Rx=vx-v1x;
        float Ry=vy-v1y;
        
        //dot products
        float RdotN=e_n[i][1]*Rx+e_n[i][2]*Ry;
        float RdotT=e_t[i][1]*Rx+e_t[i][2]*Ry;
        
        //FORCES
        if(RdotN<0 && RdotN>-std::max(ss_lA_ci, ss_lA_cj)*lALERT && RdotT>0 && RdotT<e_length[i]){
	    //printf("core overlapping");            
            float l2=pow(e_length[i],2);
            float l4=l2*l2;
            
            //v1
            v_F[v1][1]-=k0*(2*(v1y-v2y)*((v1x-v2x)*(v2x-vx)+(v1y-v2y)*(v2y-vy))*(-(v2y*vx)+v1y*(-v2x+vx)+v1x*(v2y-vy)+v2x*vy))/l4;
            v_F[v1][2]-=k0*(-2*(v1x-v2x)*((v1x-v2x)*(v2x-vx)+(v1y-v2y)*(v2y-vy))*(-(v2y*vx)+v1y*(-v2x+vx)+v1x*(v2y-vy)+v2x*vy))/l4;
            //v2
            v_F[v2][1]-=k0*(2*(v1y-v2y)*(v1y*(v2x-vx)+v2y*vx-v2x*vy+v1x*(-v2y+vy))*(pow(v1x,2)+pow(v1y,2)+v2x*vx-v1x*(v2x+vx)+v2y*vy-v1y*(v2y+vy)))/l4;
            v_F[v2][2]-=k0*(2*(v1x-v2x)*(-(v2y*vx)+v1y*(-v2x+vx)+v1x*(v2y-vy)+v2x*vy)*(pow(v1x,2)+pow(v1y,2)+v2x*vx-v1x*(v2x+vx)+v2y*vy-v1y*(v2y+vy)))/l4;
            //i
            v_F[j][1]-=k0*(2*(-v1y+v2y)*(v1y*(v2x-vx)+v2y*vx-v2x*vy+v1x*(-v2y+vy)))/l2;
            v_F[j][2]-=k0*(2*(v1x-v2x)*(v1y*(v2x-vx)+v2y*vx-v2x*vy+v1x*(-v2y+vy)))/l2;
            c_growth[cID] = 0;
	    c_growth[cell_j] = 0;
	    if (c_shrink[cID] == 0) {
	        c_A0[cID] = c_A0[cID] * 0.85;
		c_shrink[cID] = 1;
	    }
	    if (c_shrink[cell_j] == 0) {
                c_A0[cell_j] = c_A0[cell_j] / 2;
		c_shrink[cell_j] = 1;
            }		
        }
	// When a point is deep inside or the target cell has irregular shape, computing repulsion using edge might be ineffective
	if(false && RdotN > -lDEEP*std::max(ss_lA_ci, ss_lA_cj) && RdotN < 0) {
	    float dd_j_ic = torus_distance(vx, vy, c_center[cID][1], c_center[cID][2]);
	    float dd_j_jc = torus_distance(vx, vy, c_center[cell_j][1], c_center[cell_j][2]);
      	    float dx = (c_center[cell_j][1] - vx) / dd_j_jc;
            float dy = (c_center[cell_j][2] - vy) / dd_j_jc;
            v_F[j][1] += kC * dx / dd_j_ic;
	    v_F[j][2] += kC * dy / dd_j_ic;	    
	}

        
        if(RdotN>0  && RdotT>0 && RdotT<e_length[i] && (RdotN<lS*ss_ci || RdotN<lS*ss_cj)){
            //printf("v edge repul\n");
	    //printf("lS * ss_ci is %f\n", lS*ss_ci);
	    //printf("lS * ss_cj is %f\n", lS*ss_cj);
	    //printf("RdotN is %f\n", RdotN);
	    float li=e_length[i];
            float l4=pow(li,4);
            if (RdotN<lS*ss_cj && celli_A >= 0.3 * celli_CV) {
                //v1
                v_F[v1][1]-=(2*kS*(v1y - v2y)*((v1x - v2x)*(v2x - vx) + (v1y - v2y)*(v2y - vy))*(lS*ss_cj*li - v2y*vx + v1y*(-v2x + vx) + v1x*(v2y - vy) + v2x*vy))/l4;
                v_F[v1][2]-=(-2*kS*(v1x - v2x)*((v1x - v2x)*(v2x - vx) + (v1y - v2y)*(v2y - vy))*(lS*ss_cj*li - v2y*vx + v1y*(-v2x + vx) + v1x*(v2y - vy) + v2x*vy))/l4;
                //v2
                v_F[v2][1]-=(-2*kS*(v1y - v2y)*(lS*ss_cj*li - v2y*vx + v1y*(-v2x + vx) + v1x*(v2y - vy) + v2x*vy)*(pow(v1x,2) + pow(v1y,2) + v2x*vx - v1x*(v2x + vx) + v2y*vy - v1y*(v2y + vy)))/l4;
                v_F[v2][2]-=(2*kS*(v1x - v2x)*(lS*ss_cj*li - v2y*vx + v1y*(-v2x + vx) + v1x*(v2y - vy) + v2x*vy)*(pow(v1x,2) + pow(v1y,2) + v2x*vx - v1x*(v2x + vx) + v2y*vy - v1y*(v2y + vy)))/l4;
	    }
	    if (RdotN<lS*ss_ci && cellj_A >= 0.3 * cellj_CV) {
	        //i
                v_F[j][1]-=(2*kS*(v1y - v2y)*(lS*ss_ci*li - v2y*vx + v1y*(-v2x + vx) + v1x*(v2y - vy) + v2x*vy))/(pow(v1x - v2x,2) + pow(v1y - v2y,2));
                v_F[j][2]-=(2*kS*(-v1x + v2x)*(lS*ss_ci*li - v2y*vx + v1y*(-v2x + vx) + v1x*(v2y - vy) + v2x*vy))/(pow(v1x - v2x,2) + pow(v1y - v2y,2));
	    }
            edge_vertex_int[j][i]=1;
        }
    }
    
    delete []dxdydz;
    
    return 0;
    
}
//****************************************************************************
