//****************************************************************************
//****************************************************************************
//******************************ALLOCATE**************************************
//****************************************************************************
//****************************************************************************
//float **e_t;
//float **e_n;
void reset_arrays(){
    force_topology_changed = true;
    
    Nv = 0; Ne = 0; Nc = 0;
    
    
    //ids
    for(int i = 1; i<=3; i++) ids[i][0]=0;

    
    for(int i = 1; i<=array_max; i++){
        
        for(int j=1; j<=array_max; j++) edge_vertex_int[i][j]=0;
        for(int j=1; j<=array_max; j++) ve_next_time[i][j]=0;
        for(int j=1; j<=array_max; j++) vv_next_time[i][j]=0;
        
        
        //exist
        for(int j=1; j<=3; j++) exist[j][i]=0;
        
        //vertices
        for(int j = 0; j <= 2; j++) v[i][j]=0;
        for(int j = 0; j <= 2; j++) vS[i][j]=0;
        v_edges[i][0]=0;
        v_cells[i][0]=0;
        //.................
        for(int j = 0; j <= 2; j++) v_F[i][j]=0;
        //.................
        
        //edges
        for(int j = 0; j <= 2; j++) e[i][j]=0;
        e_cells[i][0]=0;
        e_Vcells[i][0]=0;
        //.................
        e_g[i]=0;
        e_length[i]=0;
        e_lengthS[i]=0;
        e_dl[i]=0;
        e_clock[i]=0;
        for(int j = 0; j <= 2; j++){
            e_t[i][j]=0;
            e_n[i][j]=0;
            e_center[i][j]=0;
        }
        
        //cells
        c[i][0]=0;
        c_vertices[i][0]=0;
        c_Vedges[i][0]=0;
        c_area[i]=0;
        c_A0[i]=0;
	c_dead[i]=0;
	c_shrink[i]=0;
        for(int j = 0; j <= 2; j++){
            c_center[i][j]=0;
        }
        c_growth[i]=0;
    }
}
void reset_c_shrink() {  
    for(int i = 1; i<=array_max; i++){
	c_shrink[i]=0;
    }
}
//****************************************************************************
void allocate(){
    
    edge_vertex_int = new int*[array_max+1];
    ve_next_time = new float*[array_max+1];
    vv_next_time = new float*[array_max+1];
    for(int i=0; i <= array_max; i++){
        edge_vertex_int[i] = new int[array_max+1];
        ve_next_time[i] = new float[array_max+1];
        vv_next_time[i] = new float[array_max+1];
    }
    
    //ids
    ids     = new int*[4];
    for(int i=0; i <= 3; i++) ids[i]   = new int[array_max+1];
    
    
    //exist
    exist     = new int*[4];
    for(int i=0; i <= 3; i++) exist[i] = new int[array_max+1];
    
    
    //vertices
    v          = new float*[array_max+1];
    vS         = new float*[array_max+1];
    v_edges   = new int*[array_max+1];
    v_cells   = new int*[array_max+1];
    //.................
    v_F        = new float*[array_max+1];
    for(int i=0; i <= array_max; i++){
        v[i]          = new float[3];
        vS[i]          = new float[3];
        v_edges[i]   = new int[200];
        v_cells[i]   = new int[200];
        //.................
        v_F[i]          = new float[3];
    }

    
    //edges
    e          = new int*[array_max+1];
    e_cells   = new int*[array_max+1];
    e_Vcells   = new int*[array_max+1];
    for(int i=0; i <= array_max; i++){
        e[i]          = new int[3];
        e_cells[i]   = new int[200];
        e_Vcells[i]   = new int[200];
    }
    e_g = new float[array_max+1];
    e_length = new float[array_max+1];
    e_lengthS = new float[array_max+1];
    e_dl = new float[array_max+1];
    e_clock = new float[array_max+1];
    e_t = new float*[array_max+1];
    e_n = new float*[array_max+1];
    e_center = new float*[array_max+1];
    for(int i=0; i <= array_max; i++){
        e_t[i] = new float[3];
        e_n[i] = new float[3];
        e_center[i] = new float[3];
    }
    
    //cells
    c   = new int*[array_max+1];
    c_dead = new int[array_max+1];
    c_shrink = new int[array_max+1];
    c_vertices   = new int*[array_max+1];
    c_Vedges   = new int*[array_max+1];
    for(int i=0; i <= array_max; i++){
        c[i]   = new int[200];
        c_vertices[i]   = new int[200];
        c_Vedges[i]   = new int[200];
    }
    c_area = new float[array_max+1];
    c_A0 = new float[array_max+1];
    c_center = new float*[array_max+1];
    for(int i=0; i <= array_max; i++){
        c_center[i] = new float[3];
    }
    c_growth = new int[array_max+1];
    
    
    
    //MISC
    perioXYZ = new float[3];
        
    
    //RESET ARRAYS
    reset_arrays();

    return ;
}
//****************************************************************************
void deallocate(){
    
    for(int i=0; i <= array_max; i++) delete [] edge_vertex_int[i];
    delete [] edge_vertex_int;
    for(int i=0; i <= array_max; i++) delete [] ve_next_time[i];
    for(int i=0; i <= array_max; i++) delete [] vv_next_time[i];
    delete [] ve_next_time;
    delete [] vv_next_time;
    
    //ids
    for(int i=0; i <= 3; i++) delete [] ids[i];
    delete [] ids;
    
    //exist
    for(int i=0; i <= 3; i++) delete [] exist[i];
    delete [] exist;
    
    
    //vertices
    for(int i=0; i <= array_max; i++){
        delete [] v[i]      ;
        delete [] vS[i]      ;
        delete [] v_edges[i]  ;
        delete [] v_cells[i]  ;
        //.................
        delete [] v_F[i]      ;
    }
    delete [] v;
    delete [] vS;
    delete [] v_edges  ;
    delete [] v_cells  ;
    //.................
    delete [] v_F;
    //.................
    
    
    //edges
    for(int i=0; i <= array_max; i++){
        delete [] e[i]      ;
        delete [] e_cells[i]  ;
        delete [] e_Vcells[i]  ;
    }
    delete [] e;
    delete [] e_cells  ;
    delete [] e_Vcells  ;
    //.................
    delete [] e_g  ;
    delete [] e_length  ;
    delete [] e_lengthS  ;
    delete [] e_dl  ;
    delete [] e_clock;
    for(int i=0; i <= array_max; i++){
        delete [] e_t[i];
        delete [] e_n[i];
        delete [] e_center[i];
    }
    delete [] e_t;
    delete [] e_n;
    delete [] e_center;
    
    //cells
    for(int i=0; i <= array_max; i++){
        delete [] c[i]  ;
        delete [] c_vertices[i]      ;
        delete [] c_Vedges[i]  ;
    }
    delete [] c  ;
    delete [] c_vertices;
    delete [] c_Vedges  ;
    delete [] c_area  ;
    delete [] c_A0  ;
    delete [] c_dead;
    delete [] c_shrink;
    for(int i=0; i <= array_max; i++){
        delete [] c_center[i];
    }
    delete [] c_center;
    delete [] c_growth  ;
    
    
    
    //MISC
    delete [] perioXYZ;

    return ;
}
//****************************************************************************
