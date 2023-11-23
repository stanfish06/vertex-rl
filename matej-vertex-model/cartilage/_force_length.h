//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
float e_length_force(int i){
            
    float ddx,ddy;
    float *dxdydz = new float[3]; dxdydz[1]=0; dxdydz[2]=0;
    int v1=e[i][1], v2=e[i][2];
    int cell_id = v_cells[i][1];
    float ss = EG[get_subdomain(cell_id)];
    torus_dx_dy_dz(dxdydz,v1,v2);
    ddx=v[v2][1]-(v[v1][1]+dxdydz[1]);
    ddy=v[v2][2]-(v[v1][2]+dxdydz[2]);
    delete []dxdydz;
    
    float c0=-e_g0*ss/e_length[i];
    
    //v1
    v_F[v1][1]+=-c0*ddx;
    v_F[v1][2]+=-c0*ddy;
    
    //v2
    v_F[v2][1]+=c0*ddx;
    v_F[v2][2]+=c0*ddy;
    
    return e_g0*ss*e_length[i];
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
