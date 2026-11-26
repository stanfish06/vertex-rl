//****************************************************************************
//****************************************************************************
//****************************T1 TRANSITION***********************************
//****************************************************************************
//****************************************************************************
int vert_edg_cells_ed(int i, int *vertices, int *edges, int *cells)
{
    //VERTICES
    int v1=e[i][1]; appendToList(vertices,v1);
    int v2=e[i][2]; appendToList(vertices,v2);
    
    //EDGES
    appendToList(edges,i);
    for(int j = 1; j<=v_edges[v1][0]; j++) appendToList(edges,v_edges[v1][j]);
    for(int j = 1; j<=v_edges[v2][0]; j++) appendToList(edges,v_edges[v2][j]);
    
    //CELLS
    int nrCells=0;
    for(int j = 1; j<=e_cells[i][0]; j++) { appendToList(cells,e_cells[i][j]); nrCells++;}
    
    return nrCells;
}
//****************************************************************************
int new_vertex(int *vertices)
{
    int vNEW=make_vertex(0,0);
    int v1=vertices[1];
    int v2=vertices[2];
    for(int k=1; k<=2; k++){
        float XYZ;
        // if(fabs(v[v2][k]-v[v1][k])>perioXYZ[k]/2.) XYZ=0.5*(v[v1][k]+v[v2][k]+perioXYZ[k]);
        // else XYZ=0.5*(v[v1][k]+v[v2][k]);
        XYZ=0.5*(v[v1][k]+v[v2][k]);
        v[vNEW][k]=XYZ;
    }
   if (v[vNEW][1] == 0) {
   	printf("Zero new vertex: \n");
	printf("Parent vertex1: (%f, %f)\n", v[v1][1], v[v1][2]);
	printf("Parent vertex1: (%f, %f)\n", v[v2][1], v[v2][2]);
   } 
    // torus_vertex(v1);
   // vertex_correction(vNEW);
    
    return vNEW;
}
//****************************************************************************
void insertEdgeIntoCell(int i, int j, int el)
{
    c[i][0]++;
    if(j==c[i][0]) c[i][j]=el;
    else {
        for(int k=c[i][0]; k>j; k--) c[i][k]=c[i][k-1];
        c[i][j]=el;
    }
}
//****************************************************************************
void correct_cells(int cID, int newe)
{
    int v1=e[newe][1]; int v2=e[newe][2];
    for(int k=1; k<=c[cID][0]; k++){
        //vprev
        int oeID=c[cID][k];
        int vprev=0;
        if(oeID>0) vprev=e[abs(oeID)][2];
        else vprev=e[abs(oeID)][1];
        //vnext
        if(k==c[cID][0]) oeID=c[cID][1];
        else oeID=c[cID][k+1];
        int vnext=0;
        if(oeID>0) vnext=e[abs(oeID)][1];
        else vnext=e[abs(oeID)][2];
        if(vprev!=vnext){
            if      (vprev==v1 && vnext==v2)   insertEdgeIntoCell(cID,k+1,newe);
            else if (vprev==v2 && vnext==v1)   insertEdgeIntoCell(cID,k+1,-newe);
            break;
        }
    }
}
//****************************************************************************
int EDGE_DIVISION(int i)
{
    
    //******************************************************************************
    int vertices[200] = {};
    int edges[200] = {};
    int cells[200] = {};
    int nrCells=vert_edg_cells_ed(i,vertices,edges,cells);
//    printf("vertices\n"); printLIST(vertices);
//    printf("edges\n"); printLIST(edges);
//    printf("cells\n"); printLIST(cells);
    //******************************************************************************


    //DISSOLVE OBJECTS
    for(int j=1; j<=cells[0]; j++) dissolve_cell(cells[j]);

    
    //CREATE NEW VERTEX
    int vNEW=new_vertex(vertices);
    
    
    //RESTITCH EDGES
    int v2PREV=e[i][2];
    remake_edge(i,v2PREV,vNEW);
    int eNEW=make_edge(vNEW,v2PREV,0);
    
    
    //CORRECT CELLS
    correct_cells(cells[1],eNEW);
    
    
    //MAKE OBJECTS
    for(int j=1; j<=cells[0]; j++) make_cell(cells[j]);

    
    //******************************************************************************
    //******************************************************************************
    
    
    return eNEW;
}
//****************************************************************************
