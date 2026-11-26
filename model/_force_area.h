//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
#include <omp.h>
float dA_new(int cellID, int j){
    
    int vert_ref_id=c_vertices[cellID][1];
    float vrefx = v[vert_ref_id][1];
    float vrefy = v[vert_ref_id][2];
    
    
    //vertices
    int v1, v2;
    if(j==c_vertices[cellID][0]){
        v1=c_vertices[cellID][j];
        v2=c_vertices[cellID][1];
    }
    else {
        v1=c_vertices[cellID][j];
        v2=c_vertices[cellID][j+1];
    }
    
    //v1
    float v1x=v[v1][1];
    float v1y=v[v1][2];
    //x
    //if(fabs(v1x-vrefx)>perioXYZ[1]/2.){
    //    if(v1x<vrefx) v1x+=perioXYZ[1];
    //    else if(v1x>vrefx) v1x-=perioXYZ[1];
    //}
    //y
    //if(fabs(v1y-vrefy)>perioXYZ[2]/2.){
    //    if(v1y<vrefy) v1y+=perioXYZ[2];
    //    else if(v1y>vrefy) v1y-=perioXYZ[2];
    //}
    
    //v2
    float v2x=v[v2][1];
    float v2y=v[v2][2];
    //x
    //if(fabs(v2x-vrefx)>perioXYZ[1]/2.){
    //    if(v2x<vrefx) v2x+=perioXYZ[1];
    //    else if(v2x>vrefx) v2x-=perioXYZ[1];
    //}
    //y
    //if(fabs(v2y-vrefy)>perioXYZ[2]/2.){
    //    if(v2y<vrefy) v2y+=perioXYZ[2];
    //    else if(v2y>vrefy) v2y-=perioXYZ[2];
    //}
    
    return 0.5*(-v2x*v1y + v1x*v2y);
}
//****************************************************************************
float CellArea_new(int i){
    float areaSum=0;
    for(int j=1; j<=c_vertices[i][0]; j++ ) areaSum += dA_new(i,j);
    return areaSum;
}
//****************************************************************************
float c_Adot_j(int i, int j){
    
    int vert_ref_id=c_vertices[i][1];
    float vrefx = v[vert_ref_id][1];
    float vrefy = v[vert_ref_id][2];
    
    //VERTEX IDs
    int v1, v2, v3;
    if(j==1){
        v1=c_vertices[i][c_vertices[i][0]];
        v2=c_vertices[i][1];
        v3=c_vertices[i][2];
    }
    else if(j==c_vertices[i][0]){
        v1=c_vertices[i][j-1];
        v2=c_vertices[i][j];
        v3=c_vertices[i][1];
    }
    else {
        v1=c_vertices[i][j-1];
        v2=c_vertices[i][j];
        v3=c_vertices[i][j+1];
    }
    
    //FORCES
    float forceX=v_F[v2][1];
    float forceY=v_F[v2][2];
    
    //VERTEX COORDINATES
    //v1
    float v1x=v[v1][1];
    float v1y=v[v1][2];
    //x
    //if(fabs(v1x-vrefx)>perioXYZ[1]/2.){
    //    if(v1x<vrefx) v1x+=perioXYZ[1];
    //    else if(v1x>vrefx) v1x-=perioXYZ[1];
    //}
    //y
    //if(fabs(v1y-vrefy)>perioXYZ[2]/2.){
    //    if(v1y<vrefy) v1y+=perioXYZ[2];
    //    else if(v1y>vrefy) v1y-=perioXYZ[2];
    //}
    
    //v3
    float v3x=v[v3][1];
    float v3y=v[v3][2];
    //x
    //if(fabs(v3x-vrefx)>perioXYZ[1]/2.){
    //    if(v3x<vrefx) v3x+=perioXYZ[1];
    //   else if(v3x>vrefx) v3x-=perioXYZ[1];
    //}
    //y
    //if(fabs(v3y-vrefy)>perioXYZ[2]/2.){
    //    if(v3y<vrefy) v3y+=perioXYZ[2];
    //    else if(v3y>vrefy) v3y-=perioXYZ[2];
    //}
    
    return forceX*(v3y-v1y)-forceY*(v3x-v1x);
}
//****************************************************************************
float c_AdotCalc(int i){
    float sum=0;
    for(int j=1; j<=c_vertices[i][0]; j++) sum+=c_Adot_j(i,j);
    return 0.5*sum;
}
//****************************************************************************
void assemble_area_force(int cell, ForceWorkspace& work) {
    int n = c_vertices[cell][0];
    float derivative = -2.*kA*(c_area[cell]-c_A0[cell]);
    for (int j = 1; j <= n; ++j) {
        int a = c_vertices[cell][j];
        int b = c_vertices[cell][j == n ? 1 : j + 1];
        // Each polygon owns its vertices. At vertex 1 the outgoing edge
        // comes first; all other vertices see the incoming edge first.
        work.area_terms[a][j == 1 ? 0 : 1] = {derivative*(v[b][2]/2.), derivative*(-v[b][1]/2.)};
        work.area_terms[b][j == n ? 1 : 0] = {derivative*(-v[a][2]/2.), derivative*(v[a][1]/2.)};
    }
    work.area_energy[cell] = kA*pow(c_area[cell]-c_A0[cell], 2);
}
