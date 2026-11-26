//****************************************************************************
//****************************************************************************
//****************************************************************************
//****************************************************************************
void edge_center(int i){
    
    int v1=e[i][1], v2=e[i][2];
    
    for(int k=1; k<=2; k++){
        float XYZ;
        // if(fabs(v[v2][k]-v[v1][k])>perioXYZ[k]/2.) XYZ=0.5*(v[v1][k]+v[v2][k]+perioXYZ[k]);
        // else XYZ=0.5*(v[v1][k]+v[v2][k]);
        XYZ=0.5*(v[v1][k]+v[v2][k]);
        e_center[i][k]=XYZ;
    }
    //TORUS
    //x
    // if      (e_center[i][1]<0)             e_center[i][1] += perioXYZ[1];
    // else if (e_center[i][1]>perioXYZ[1])   e_center[i][1] -= perioXYZ[1];
    //y
    // if      (e_center[i][2]<0)             e_center[i][2] += perioXYZ[2];
    // else if (e_center[i][2]>perioXYZ[2])   e_center[i][2] -= perioXYZ[2];
}
//****************************************************************************
int e_sign_cell(int i, int cID){
    int sgn=0;
    for(int j=1; j<=c[cID][0]; j++){
        int oeID=c[cID][j];
        int eID=abs(oeID);
        if(eID==i){
            if(oeID<0) sgn=-1;
            else sgn=1;
            break;
        }
    }
    return sgn;
}
//****************************************************************************
void tangent_and_normal(int i){
    
    int v1=e[i][1], v2=e[i][2];
    int cID=e_cells[i][1];
    int sgn=e_sign_cell(i,cID);
    if(sgn<0){
        int vv=v1;
        v1=v2;
        v2=vv;
    }
    
    float *dxdydz = new float[3]; dxdydz[1]=0; dxdydz[2]=0;
    torus_dx_dy_dz(dxdydz,v1,v2);
    float ddx=v[v2][1]-(v[v1][1]+dxdydz[1]);
    float ddy=v[v2][2]-(v[v1][2]+dxdydz[2]);
    delete []dxdydz;
    float elen=e_length[i];
    
    e_t[i][1]=ddx/elen;
    e_t[i][2]=ddy/elen;
    e_n[i][1]=e_t[i][2];
    e_n[i][2]=-e_t[i][1];
}
//****************************************************************************
float edge_length(int i){
    int v1=e[i][1], v2=e[i][2];
    float *dxdydz = new float[3]; dxdydz[1]=0; dxdydz[2]=0;
    torus_dx_dy_dz(dxdydz,v1,v2);
    float ddx=v[v2][1]-(v[v1][1]+dxdydz[1]);
    float ddy=v[v2][2]-(v[v1][2]+dxdydz[2]);
    delete []dxdydz;
    return sqrt(ddx*ddx+ddy*ddy);
}
//****************************************************************************
float vv_distance(int i, int j){
    float *dxdydz = new float[3]; dxdydz[1]=0; dxdydz[2]=0;
    torus_dx_dy_dz(dxdydz,i,j);
    float ddx=v[j][1]-(v[i][1]+dxdydz[1]);
    float ddy=v[j][2]-(v[i][2]+dxdydz[2]);
    delete []dxdydz;
    return sqrt(ddx*ddx+ddy*ddy);
}
//****************************************************************************
