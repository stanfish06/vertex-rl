//****************************************************************************
//****************************************************************************
//*******************************TORUS****************************************
//****************************************************************************
//****************************************************************************
void torus_vertex(int i){
    //x
    //if      (v[i][1]<0)             v[i][1] += perioXYZ[1];
    //else if (v[i][1]>perioXYZ[1])   v[i][1] -= perioXYZ[1];
    //y
    //if      (v[i][2]<0)             v[i][2] += perioXYZ[2];
    //else if (v[i][2]>perioXYZ[2])   v[i][2] -= perioXYZ[2];

    if      (v[i][1]<0)             v[i][1] = 0;
    else if (v[i][1]>perioXYZ[1])   v[i][1] = perioXYZ[1];
    //y
    if      (v[i][2]<0)             v[i][2] = 0;
    else if (v[i][2]>perioXYZ[2])   v[i][2] = perioXYZ[2];
}
//****************************************************************************
void torus_dx_dy_dz(float *dxdydz, int vert_id, int vert_ref_id){

    dxdydz[1]=0; dxdydz[2]=0;

    //x
    //if(fabs(v[vert_id][1]-v[vert_ref_id][1])>0.5*perioXYZ[1]){
    //    if      (v[vert_id][1]<v[vert_ref_id][1])   dxdydz[1]=perioXYZ[1];
    //    else if (v[vert_id][1]>v[vert_ref_id][1])   dxdydz[1]=-perioXYZ[1];
    //}
    //else dxdydz[1]=0;
    //y
    //if(fabs(v[vert_id][2]-v[vert_ref_id][2])>0.5*perioXYZ[2]){
    //    if      (v[vert_id][2]<v[vert_ref_id][2])   dxdydz[2]=perioXYZ[2];
    //    else if (v[vert_id][2]>v[vert_ref_id][2])   dxdydz[2]=-perioXYZ[2];
    //}
    //else dxdydz[2]=0;
}
//****************************************************************************
void torus_dx_dy_dzS(float *dxdydzS, int vert_id, int vert_ref_id){

    dxdydzS[1]=0; dxdydzS[2]=0;

    //x
    if(fabs(vS[vert_id][1]-vS[vert_ref_id][1])>0.5*perioXYZ[1]){
        if      (vS[vert_id][1]<vS[vert_ref_id][1])   dxdydzS[1]=perioXYZ[1];
        else if (vS[vert_id][1]>vS[vert_ref_id][1])   dxdydzS[1]=-perioXYZ[1];
    }
    else dxdydzS[1]=0;
    //y
    if(fabs(vS[vert_id][2]-vS[vert_ref_id][2])>0.5*perioXYZ[2]){
        if      (vS[vert_id][2]<vS[vert_ref_id][2])   dxdydzS[2]=perioXYZ[2];
        else if (vS[vert_id][2]>vS[vert_ref_id][2])   dxdydzS[2]=-perioXYZ[2];
    }
    else dxdydzS[2]=0;
}
//****************************************************************************
void torus_dx_dy_dzS2(float *dxdydzS, int vert_id, int vert_ref_id){

    dxdydzS[1]=0; dxdydzS[2]=0;

    //x
    if(fabs(vS[vert_id][1]-v[vert_ref_id][1])>0.5*perioXYZ[1]){
        if      (vS[vert_id][1]<v[vert_ref_id][1])   dxdydzS[1]=perioXYZ[1];
        else if (vS[vert_id][1]>v[vert_ref_id][1])   dxdydzS[1]=-perioXYZ[1];
    }
    else dxdydzS[1]=0;
    //y
    if(fabs(vS[vert_id][2]-v[vert_ref_id][2])>0.5*perioXYZ[2]){
        if      (vS[vert_id][2]<v[vert_ref_id][2])   dxdydzS[2]=perioXYZ[2];
        else if (vS[vert_id][2]>v[vert_ref_id][2])   dxdydzS[2]=-perioXYZ[2];
    }
    else dxdydzS[2]=0;
}
//****************************************************************************
float minc(float a1, float a2){
    if(a1<a2) return a1;
    else return a2;
}
//****************************************************************************
float maxc(float a1, float a2){
    if(a1<a2) return a2;
    else return a1;
}
//****************************************************************************
float torus_distance(float xc, float yc, float xr, float yr){
    //x
    float a1=fabs(xc-xr);
   // float mi=minc(xc,xr);
   // float ma=maxc(xc,xr);
   // float a2=fabs(mi+perioXYZ[1]-ma);
   // float rx=minc(a1,a2);
    //y
    float a2=fabs(yc-yr);
   // mi=minc(yc,yr);
   // ma=maxc(yc,yr);
   // a2=fabs(mi+perioXYZ[2]-ma);
   // float ry=minc(a1,a2);
    return sqrt(a1*a1+a2*a2);
}
//****************************************************************************
