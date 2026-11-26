//****************************************************************************
//****************************************************************************
//******************************OUTPUT****************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
float out_polygons(int fcount){
    
    
    char filename2[100]; char filename3[100];
    snprintf(filename2, sizeof(char) * 100, "./output/X_%d.txt", fcount);
    snprintf(filename3, sizeof(char) * 100, "./output/Y_%d.txt", fcount);
    
    float *dxdydz = new float[3];
    dxdydz[1]=0; dxdydz[2]=0;
    int vert_ref_id, vert_id;
    FILE *file2; file2 = fopen(filename2, "wt");
    FILE *file3; file3 = fopen(filename3, "wt");
    
    for(int i=1; i<=Nc; i++) if(exist[3][i]==1){
        
        int srgn=1;
            
        fprintf(file2, "%d ", srgn*c_vertices[i][0]);
        fprintf(file3, "%d ", srgn*c_vertices[i][0]);
            
        vert_ref_id=c_vertices[i][1];
        fprintf(file2, "%.16g ", v[vert_ref_id][1]);
        fprintf(file3, "%.16g ", v[vert_ref_id][2]);
            
        for(int j = 2; j <= c_vertices[i][0]; j++){
            vert_id=c_vertices[i][j];
            torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
            fprintf(file2, "%.16g ", v[vert_id][1] + dxdydz[1]);
            fprintf(file3, "%.16g ", v[vert_id][2] + dxdydz[2]);
        }
            
        for(int j = c_vertices[i][0]; j <= 15; j++){
            fprintf(file2, "0 ");
            fprintf(file3, "0 ");
        }
            
        fprintf(file2, "\n");
        fprintf(file3, "\n");
    }
    fclose(file2);
    fclose(file3);
    
    //SHOULDER
    char filename4[100]; char filename5[100];
    snprintf(filename4, sizeof(char) * 100, "./output/XS_%d.txt", fcount);
    snprintf(filename5, sizeof(char) * 100, "./output/YS_%d.txt", fcount);
    FILE *file4; file4 = fopen(filename4, "wt");
    FILE *file5; file5 = fopen(filename5, "wt");
    for(int i=1; i<=Nc; i++) if(exist[3][i]==1){
        
        int srgn=1;
            
        fprintf(file4, "%d ", srgn*c_vertices[i][0]);
        fprintf(file5, "%d ", srgn*c_vertices[i][0]);
            
        vert_ref_id=c_vertices[i][1];
        fprintf(file4, "%.16g ", vS[vert_ref_id][1]);
        fprintf(file5, "%.16g ", vS[vert_ref_id][2]);
            
        for(int j = 2; j <= c_vertices[i][0]; j++){
            vert_id=c_vertices[i][j];
            //torus_dx_dy_dzS(dxdydz,vert_id,vert_ref_id);
            fprintf(file4, "%.16g ", vS[vert_id][1]);
            fprintf(file5, "%.16g ", vS[vert_id][2]);
        }
            
        for(int j = c_vertices[i][0]; j <= 15; j++){
            fprintf(file4, "0 ");
            fprintf(file5, "0 ");
        }
            
        fprintf(file4, "\n");
        fprintf(file5, "\n");
    }
    fclose(file4);
    fclose(file5);
    
    delete []dxdydz;
    
    return 0;
}

float out_to_df(int fcount) {
    char filename[100]; char filename2[100];
    snprintf(filename, sizeof(char) * 100, "./output/snapshot_%d.txt", fcount);
    //snprintf(filename2, sizeof(char) * 100, "./initial/last_snapshot_%d.txt", fcount);
    FILE *file; file = fopen(filename, "wt");
    //FILE *file2; file2 = fopen(filename2, "wt");
    fprintf(file, "cell_id,v_x,v_y,vs_x,vs_y,c_x,c_y,A,A0,dead\n");
    //fprintf(file2, "cell_id,v_x,v_y,vs_x,vs_y,c_x,c_y,A,A0,dead\n");
    for (int i = 1; i <= Nc; i++) if(exist[3][i] == 1){
	
	float c_x = c_center[i][1];
	float c_y = c_center[i][2];
	float cell_A = c_area[i];
	float cell_A0 = c_A0[i];
	int cell_dead = c_dead[i];
	
        for(int j = 1; j <= c_vertices[i][0]; j++){
            int vert_id=c_vertices[i][j];
	    float v_x = v[vert_id][1];
	    float v_y = v[vert_id][2];
	    // Shoulder v xy is flipped for some reason, fix it later
	    float vs_x = vS[vert_id][2];
	    float vs_y = vS[vert_id][1];
            fprintf(file, "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%d\n", i, v_x, v_y, vs_x, vs_y, c_x, c_y, cell_A, cell_A0, cell_dead);
           // fprintf(file2, "%d,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%.16g,%d\n", i, v_x, v_y, vs_x, vs_y, c_x, c_y, cell_A, cell_A0, cell_dead);
       }
       
    }
    fclose(file);
    //fclose(file2);
    return 0;
}


float out_cell_centroid(int fcount) {
    char filename2[100]; char filename3[100];
    snprintf(filename2, sizeof(char) * 100, "./output/XC_%d.txt", fcount);
    snprintf(filename3, sizeof(char) * 100, "./output/YC_%d.txt", fcount);
    
    //float *dxdydz = new float[3];
    //dxdydz[1]=0; dxdydz[2]=0;
    int vert_ref_id, vert_id;
    FILE *file2; file2 = fopen(filename2, "wt");
    FILE *file3; file3 = fopen(filename3, "wt");
    for(int i=1; i<=Nc; i++) if(exist[3][i]==1){
        fprintf(file2, "%d ", i);
        fprintf(file3, "%d ", i);

        fprintf(file2, "%.16g ", c_center[i][1]);
        fprintf(file3, "%.16g ", c_center[i][2]);
            
        fprintf(file2, "\n");
        fprintf(file3, "\n");
    }
    fclose(file2);
    fclose(file3);

    return 0;
}


float out_cell_state(int fcount) {
    char filename2[100];
    snprintf(filename2, sizeof(char) * 100, "./output/C_dead_%d.txt", fcount);
    
    FILE *file2; file2 = fopen(filename2, "wt");
    for(int i=1; i<=Nc; i++) if(exist[3][i]==1){
        fprintf(file2, "%d ", i);

        fprintf(file2, "%d ", c_dead[i]);
            
        fprintf(file2, "\n");
    }
    fclose(file2);

    return 0;
}





//****************************************************************************
//float out_edges(int fcount){
//
//    char filename2[50];
//    snprintf(filename2, sizeof(char) * 50, "./output/edges_%g_%d_%d_%d_%d.txt", P0, oneOvertauINT, RhoINT, THEedge, fcount);
//
//    float *dxdydz = new float[2]; dxdydz[0]=0; dxdydz[1]=0;
//    FILE *file2; file2 = fopen(filename2, "wt");
//
//    for(int i=1; i<=Ne; i++) if(e[i][0]!=0){
//        int vert_ref_id=e[i][1];
//        int vert_id=e[i][2];
//        torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
//        fprintf(file2, "%g  %g  %g  %g\n", v[vert_ref_id][1], v[vert_ref_id][2], v[vert_id][1] + dxdydz[0], v[vert_id][2] + dxdydz[1]);
//    }
//    fclose(file2);
//
//    delete []dxdydz;
//
//    return 0;
//}
//****************************************************************************
void out_tissue(int fcount){
    //out_edges(fcount);
    out_polygons(fcount);
    out_cell_centroid(fcount);
    out_cell_state(fcount);
    out_to_df(fcount);
}
//****************************************************************************
//void out_Vertissue2D(){
//    
//    char filename2[100];
//    snprintf(filename2, sizeof(char) * 100, "./output/final_%g.vt2d", P0);
//    FILE *file1; file1 = fopen(filename2, "wt");
//    
//    int nrV=0, nrE=0, nrC=0;
//    
//    //ALLOCATE INTERNAL VARIABLES
//    int **v_natural_id; int **e_natural_id; int **c_natural_id;
//    v_natural_id  = new int*[Nv+1]; e_natural_id  = new int*[Ne+1]; c_natural_id  = new int*[Nc+1];
//    for(int i=0; i < Nv+1; i++) v_natural_id[i]  = new int[3];
//    for(int i=0; i < Ne+1; i++) e_natural_id[i]  = new int[3];
//    for(int i=0; i < Nc+1; i++) c_natural_id[i]  = new int[3];
//    
//    //RESET NATURAL IDS
//    for(int i=1; i<=Nv; i++){
//        v_natural_id[i][1]=0;
//        v_natural_id[i][2]=0;
//    }
//    for(int i=1; i<=Ne; i++){
//        e_natural_id[i][1]=0;
//        e_natural_id[i][2]=0;
//    }
//    for(int i=1; i<=Nc; i++){
//        c_natural_id[i][1]=0;
//        c_natural_id[i][2]=0;
//    }
//    
//    //VERTICES
//    int cnt=0;
//    for(int i=1; i<=Nv; i++) if(exist[1][i]==1){
//        cnt++;
//        v_natural_id[i][1]=cnt;
//        v_natural_id[cnt][2]=i;
//    }
//    nrV=cnt;
//    
//    //EDGES
//    cnt=0;
//    for(int i=1; i<=Ne; i++) if(exist[2][i]==1){
//        cnt++;
//        e_natural_id[i][1]=cnt;
//        e_natural_id[cnt][2]=i;
//    }
//    nrE=cnt;
//    for(int i=1; i<=Nv; i++) if(exist[1][i]==1) fprintf(file1, "%f  %f\n", v[i][1], v[i][2]);
//    fprintf(file1, "\n");
//    
//    //EDGES
//    for(int i=1; i<=Ne; i++) if(exist[2][i]==1) fprintf(file1, "%d  %d\n", v_natural_id[e[i][1]][1], v_natural_id[e[i][2]][1]);
//    fprintf(file1, "\n");
//    
//    //CELLS
//    for(int i=1; i<=Nc; i++) if(exist[3][i]==1){
//        fprintf(file1, "%d  ", c[i][0]);
//        for(int j=1; j<=c[i][0]; j++){
//            if(c[i][j]<0) fprintf(file1, "%d  ", -e_natural_id[abs(c[i][j])][1]);
//            else fprintf(file1, "%d  ", e_natural_id[abs(c[i][j])][1]);
//        }
//        for(int j=15; j>c[i][0]; j--){
//            fprintf(file1, "0  ");
//        }
//        fprintf(file1, "\n");
//    }
//    
//    //DEALLOCATE INTERNAL VARIABLES
//    for(int i=0; i < Nv+1; ++i) delete [] v_natural_id[i] ;
//    for(int i=0; i < Ne+1; ++i) delete [] e_natural_id[i] ;
//    for(int i=0; i < Nc+1; ++i) delete [] c_natural_id[i] ;
//    delete [] v_natural_id;
//    delete [] e_natural_id;
//    delete [] c_natural_id;
//    
//    fclose(file1);
//}
//****************************************************************************
