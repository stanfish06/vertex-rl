//****************************************************************************
//****************************************************************************
//****************************T1 TRANSITION***********************************
//****************************************************************************
//****************************************************************************
int vert_edg_cells_em(int i, int *vertices, int *edges, int *cells)
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
void merge_vertices(int *vertices)
{
    int v1=vertices[1];
    int v2=vertices[2];
    for(int k=1; k<=2; k++){
        float XYZ;
        // if(fabs(v[v2][k]-v[v1][k])>perioXYZ[k]/2.) XYZ=0.5*(v[v1][k]+v[v2][k]+perioXYZ[k]);
        // else XYZ=0.5*(v[v1][k]+v[v2][k]);
        XYZ=0.5*(v[v1][k]+v[v2][k]);
        v[v1][k]=XYZ;
    }
    // torus_vertex(v1);
}
//****************************************************************************
int EDGE_MERGE(int i)
{
    
    //******************************************************************************
    int *vertices; vertices = new int[200]; for(int j=0; j<=199; j++) vertices[j]=0;
    int *edges; edges = new int[200]; for(int j=0; j<=199; j++) edges[j]=0;
    int *cells; cells = new int[200]; for(int j=0; j<=199; j++) cells[j]=0;
    int nrCells=vert_edg_cells_em(i,vertices,edges,cells);
    //******************************************************************************


    //DISSOLVE OBJECTS
    for(int j=1; j<=cells[0]; j++) dissolve_cell(cells[j]);


    //REMOVE EDGE i FROM LIST OF CELLS
    for(int j=1; j<=nrCells; j++) removeELEMENT(c,cells[j],i);


    //DISSOLVE EDGE
    dissolve_edge(i,0);


    //MERGE VERTEX POSITIONS
    merge_vertices(vertices);


    //RESTITCH EDGES
    remake_edge(edges[edges[0]],vertices[2],vertices[1]);


    //DISSOLVE VERTEX v2
    dissolve_vertex(vertices[2]);


    //MAKE OBJECTS
    for(int j=1; j<=cells[0]; j++) make_cell(cells[j]);

    
    //******************************************************************************
    delete [] vertices;
    delete [] edges;
    delete [] cells;
    //******************************************************************************
    
    
    return 0;
}
//****************************************************************************
