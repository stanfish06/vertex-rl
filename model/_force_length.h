//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void assemble_length_force(int i, ForceWorkspace& work) {
    int a = e[i][1], b = e[i][2];
    // Dev uses the vertex slot with this edge's index for the tension zone.
    float ss = EG[get_subdomain(v_cells[i][1])];
    float c0 = -e_g0*ss/e_length[i];
    float dx = v[b][1] - v[a][1], dy = v[b][2] - v[a][2];
    work.length_terms[i][0] = {-c0*dx, -c0*dy};
    work.length_terms[i][1] = {c0*dx, c0*dy};
    work.length_energy[i] = e_g0*ss*e_length[i];
}
//****************************************************************************
float e_ldotCalc(int i){
    
    float ddx,ddy;
    float *dxdydz = new float[3]; dxdydz[1]=0; dxdydz[2]=0;
    int v1=e[i][1], v2=e[i][2];
    torus_dx_dy_dz(dxdydz,v1,v2);
    ddx=v[v2][1]-(v[v1][1]+dxdydz[1]);
    ddy=v[v2][2]-(v[v1][2]+dxdydz[2]);
    delete []dxdydz;
    
    return (ddx*(v_F[v2][1]-v_F[v1][1])+ddy*(v_F[v2][2]-v_F[v1][2]))/e_length[i];
}
//****************************************************************************
