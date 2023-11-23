//****************************************************************************
//****************************************************************************
//******************************DISSOLVE**************************************
//****************************************************************************
//****************************************************************************
void dissolve_vertex(int i){

    //ID
    freeID(1,i);

    //DISSOLVE
    exist[1][i]=0;
    /* v_cells */   v_cells[i][0]=0;
    /* v_edges */   v_edges[i][0]=0;
    /* v */         v[i][0]=0;
    
    //ATTRIBUTES
    for(int j=1; j<=2; j++) v_F[i][j]=0;
    for(int j=1; j<=Ne; j++) if(exist[2][j]==1) ve_next_time[j][i]=0;
    for(int j=1; j<=Nv; j++) if(exist[1][j]==1){
        vv_next_time[i][j]=0;
        vv_next_time[j][i]=0;
    }
    
}
//****************************************************************************
void reset_edge_attributes(int i){
    e_g[i]=0;
    e_length[i]=0;
    e_dl[i]=0;
    e_clock[i]=-1;
    for(int j=1; j<=Nv; j++) if(exist[1][j]==1) ve_next_time[i][j]=0;
}
//****************************************************************************
int dissolve_edge(int i, int reserved){
    //ID
    if(reserved==0) freeID(2,i);

    //EXTERNAL ATTRIBUTES
    /* v_edges */   for(int j=1; j<=2; j++) {
	//printf("start removing edge: %d\n", i); 
	removeELEMENT(v_edges,e[i][j],i);
    }
    //DISSOLVE
    exist[2][i]=0;
    /* e_Vcells */  e_Vcells[i][0]=0;
    /* e_cells */   e_cells[i][0]=0;
    /* e */         e[i][0]=0;
    
    //ATTRIBUTES
    if(reserved==0) reset_edge_attributes(i);
    
    return i;
}
//****************************************************************************
void dissolve_cell(int i){
    //EXTERNAL ATTRIBUTES
    /* e_Vcells*/   for(int j=1; j<=c_Vedges[i][0]; j++)    {/*printf("start removing e_Vcells: %d\n", i);*/  removeELEMENT(e_Vcells,c_Vedges[i][j],i);}
    /* e_cells */   for(int j=1; j<=c[i][0]; j++)           {/*printf("start removing e_cells: %d\n", i);*/  removeELEMENT(e_cells,abs(c[i][j]),i);}
    /* v_cells */   for(int j=1; j<=c_vertices[i][0]; j++)  {/*printf("start removing v_cells: %d\n", i);*/  removeELEMENT(v_cells,c_vertices[i][j],i);}

    //DISSOLVE
    exist[3][i]=0;
    /* c_Vedges */  c_Vedges[i][0]=0;
    /* c_vertices */c_vertices[i][0]=0;

}
//****************************************************************************
void delete_cell(int i){

    //ID
    freeID(3,i);

    //DISSOLVE
    c[i][0]=0;
    c_area[i]=0;
    c_A0[i]=0;
}
//****************************************************************************
