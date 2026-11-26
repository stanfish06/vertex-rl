//****************************************************************************
//****************************************************************************
//****************************T1 TRANSITION***********************************
//****************************************************************************
//****************************************************************************
#include <omp.h>
void calc_cell_center(int i){
    //float *dxdydz = new float[3];
    //dxdydz[1]=0; dxdydz[2]=0;
    c_center[i][1]=0; c_center[i][2]=0;
    //ref vertex
    //int vert_ref_id=c_vertices[i][1];
    //c_center[i][1]+=v[vert_ref_id][1];
    //c_center[i][2]+=v[vert_ref_id][2];
    //other vertices
    for(int j=1; j<=c_vertices[i][0]; j++){
        int vert_id=c_vertices[i][j];
        //torus_dx_dy_dz(dxdydz,vert_id,vert_ref_id);
        c_center[i][1]+=v[vert_id][1];
        c_center[i][2]+=v[vert_id][2];
    }
    c_center[i][1]/=1.*c_vertices[i][0];
    c_center[i][2]/=1.*c_vertices[i][0];
    //delete []dxdydz;
}
//****************************************************************************
// Opt-in qualification of the episode RNG stream. Only the serial fate pass
// calls this function; the timed runs leave VERTEX_FATE_TRACE unset.
void trace_fate(int cell, const char* event, float draw, float outcome) {
    static const char* path = getenv("VERTEX_FATE_TRACE");
    if (!path) return;
    static std::unique_ptr<FILE, decltype(&fclose)> file([&]() {
        FILE* stream = fopen(path, "w");
        if (!stream) throw std::runtime_error("cannot open fate trace");
        fprintf(stream, "time,cell_id,event,draw,outcome\n");
        return stream;
    }(), fclose);
    if (fprintf(file.get(), "%.9g,%d,%s,%.9g,%.9g\n", Time, cell, event, draw, outcome) < 0) {
        throw std::runtime_error("cannot write fate trace");
    }
}

bool fate_draw(int cell, const char* event, double threshold) {
    float draw = rnd();
    bool accepted = threshold > draw;
    trace_fate(cell, event, draw, accepted);
    return accepted;
}

int DIVIDE_CELL(int i, float separ, int subdomain)
{
    //printf("Start dividing");
    calc_cell_center(i);
    
    float ddx,ddy;
    
    float kk;
    //float angl=M_PI*rnd();
    //float angl=M_PI;
    float c_y = c_center[i][2];
    if (subdomain == 1) printf("Proliferative zone division\n");
    float draw = rnd();
    float angl = M_PI * draw;
    trace_fate(i, "angle", draw, angl);
    //printf("cell: %d", i);
    printf("angle: %f\n", angl);
    if(angl>M_PI/2.){
        angl=M_PI-angl;
        kk=-tan(angl);
    }
    else kk=tan(angl);
    
    float cx=c_center[i][1];
    float cy=c_center[i][2];
    
    float nn=-kk*cx+cy;
    
    // c[i][0] is the number of sides of each cell
    int j1=0, j2=0;
    for(int j=1; j<=c[i][0]; j++){
        int eID=abs(c[i][j]);
        int v1=e[eID][1], v2=e[eID][2];
        //torus_dx_dy_dz(dxdydz,v1,v2);
        float v1x=v[v1][1];
        float v1y=v[v1][2];
        float v2x=v[v2][1];
        float v2y=v[v2][2];
        float dy1=kk*v1x+nn-v1y;
        float dy2=kk*v2x+nn-v2y;
        if(dy1*dy2<0){
            if(j1==0) j1=j;
            else j2=j;
        }
    }
    
    if (!j1 || !j2) throw std::runtime_error("division line misses the polygon");
    int jj;
    if(j1>j2){
        jj=j1;
        j1=j2;
        j2=jj;
    }
    
    //******************************************************************************
    int c_copy[200] = {};
    //******************************************************************************

    // This step is figuring out edge for division
    int e1=abs(c[i][j1]); int v11=e[e1][1]; int v12=e[e1][2];
    int e2=abs(c[i][j2]); int v21=e[e2][1]; int v22=e[e2][2];

    float ec1x=e_center[e1][1];
    float ec1y=e_center[e1][2];
    float ec2x=e_center[e2][1];
    float ec2y=e_center[e2][2];

    //DISSOLVE OBJECTS
    dissolve_cell(i);

    //MAKE NEW VERTICES AND EDGES
    //   INITIAL
    //   (v12)o<-----------------(e1)-------------o(v11)
    //   (v21)o------------------(e2)------------>o(v22)

    //   FINAL
    //   (v12)o<--(e1NEW)--o(ve12)(ve11)o<--(e1)--o(v11)
    //                     |            |
    //                     |            |
    //                     |            |
    //                    (ed2)        (ed1)
    //                     |            |
    //                     |            |
    //                     |            |
    //                     v            v
    //   (v21)o---(e2)---->o(ve21)(ve22)o---(e2NEW)-->o(v22)

    //e1
    //torus_dx_dy_dz(dxdydz,v11,v12);
    ddx=v[v12][1]-(v[v11][1]); ddx/=e_length[e1];
    ddy=v[v12][2]-(v[v11][2]); ddy/=e_length[e1];
    int ve11=make_vertex(ec1x-separ*ddx,ec1y-separ*ddy);
    int ve12=make_vertex(ec1x+separ*ddx,ec1y+separ*ddy);
    //e2
    //torus_dx_dy_dz(dxdydz,v21,v22);
    ddx=v[v22][1]-(v[v21][1]); ddx/=e_length[e2];
    ddy=v[v22][2]-(v[v21][2]); ddy/=e_length[e2];
    int ve21=make_vertex(ec2x-separ*ddx,ec2y-separ*ddy);
    int ve22=make_vertex(ec2x+separ*ddx,ec2y+separ*ddy);

    int ed1, ed2;
    if(c[i][j1]>0 && c[i][j2]>0){
        ed1=make_edge(ve11,ve22,0);
        ed2=make_edge(ve12,ve21,0);
    }
    if(c[i][j1]<0 && c[i][j2]<0){
        ed1=make_edge(ve12,ve21,0);
        ed2=make_edge(ve11,ve22,0);
    }
    if(c[i][j1]>0 && c[i][j2]<0){
        ed1=make_edge(ve11,ve21,0);
        ed2=make_edge(ve12,ve22,0);
    }
    if(c[i][j1]<0 && c[i][j2]>0){
        ed1=make_edge(ve12,ve22,0);
        ed2=make_edge(ve11,ve21,0);
    }


    //CORRECT EDGES
    remake_edge(e1,v12,ve11);
    int e1NEW=make_edge(ve12,v12,0);
    remake_edge(e2,v22,ve21);
    int e2NEW=make_edge(ve22,v22,0);

    //cell 1
    if(c[i][j1]>0) appendToList(c_copy,e1NEW);
    else appendToList(c_copy,c[i][j1]);
    for(int j=j1+1; j<j2; j++) appendToList(c_copy,c[i][j]);
    if(c[i][j2]>0) appendToList(c_copy,c[i][j2]);
    else appendToList(c_copy,-e2NEW);
    appendToList(c_copy,-ed2);

    //cell 2
    int sides[200] = {};
    int j=j2;
    if(c[i][j2]>0) appendToList(sides,e2NEW);
    else appendToList(sides,c[i][j2]);
    while(1<2){
        j++;
        if(j>c[i][0]) j=1;
        if(j==j1) break;
        else appendToList(sides,c[i][j]);
    }
    if(c[i][j1]>0) appendToList(sides,c[i][j1]);
    else appendToList(sides,-e1NEW);
    appendToList(sides,ed1);


    for(int j=0; j<=c_copy[0]; j++) c[i][j]=c_copy[j];

    int cIDNEW=cell_sides(sides[0],sides,0);

    //REMAKE CELLS
    make_cell(i);
    make_cell(cIDNEW);


    //******************************************************************************
    //******************************************************************************

    
    return cIDNEW;
}
//****************************************************************************

void growth_division(){
    float nrSUSC=0;
//#pragma omp parallel for num_threads(num_threads_count) reduction(+: nrSUSC) 
    //for(int i=1; i<=Nc; i++) if(exist[3][i]==1) if(get_subdomain(i) == 1 && c_area[i] > CV[get_subdomain(i)]) {
    //    nrSUSC++;
    //}
    //nrSUSC = nrSUSC / proliferative_zone_area;
   // printf("Number of dividing cells: %f\n", nrSUSC);
    // New daughters enter the fate pass on the following substep.
    const int initial_cells = Nc;
    for(int i=1; i<=initial_cells; i++) if(exist[3][i]==1){
        //GROWTH
        //if(c_growth[i]==0 && rand() < 0.5) c_growth[i]=1;
        //if(c_growth[i]==0 && rnd()<0.00001/(nrSUSC*1.)) c_growth[i]=1;
	//DIVISION
	// cell_volume is the volume needed for division
	int cell_domain = get_subdomain(i);
        float cell_volume = CV[cell_domain];
	//float c_y = c_center[i][2];
        //if (c_y >= resting_zone_interface[0] && c_y <= resting_zone_interface[1]) {
        //    cell_volume = 1;
        //} else if (c_y >= hypertrophic_zone_interface[0] && c_y <= hypertrophic_zone_interface[1]) {
        //    cell_volume = 4;
        //} else if (c_y >= proliferative_zone_interface[0] && c_y <= proliferative_zone_interface[1]) {
        //    cell_volume = 2;
        //} else {
	//    cell_volume = 1;
	//}
	// maximum desired volume is slightly bigger than the minimum division volume as optimal area is going to be slightly smaller than the desired area due to other forces
        if(c_growth[i]==0 && c_dead[i] == 0 && c_A0[i] <= cell_volume * 1.5) c_growth[i]=1;
        else if (c_growth[i]==1 && c_A0[i] <= cell_volume * 1.5) c_A0[i]+=h * 5 * (1 - c_A0[i]/(1.5*cell_volume));
//	else if (c_growth[i]==1 && c_A0[i] > cell_volume * 1.5) c_growth[i] = 0;
        //printf("Cell area: %f", c_A0[i]);
	if (cell_domain == 3 && c_area[i] >= 0.8 * cell_volume && c_dead[i] == 0) {
	    if (fate_draw(i, "death", h * 1/10)) {
	        c_dead[i] = 1;
		c_growth[i] = 0;
		printf("Cell %d just died!\n", i);
	    }
	}

	// desired area will increase if it is smaller than 1.5 fold of the maximum area
        if(c_area[i]>cell_volume && c[i][0]>=min_division_n_edges) {
	    if (cell_domain == 1 && fate_draw(i, "divide", h * 0.8)) {
               printf("start division");
               int cNEW=DIVIDE_CELL(i,0.001,cell_domain);
              // c_A0[i]=1; c_A0[cNEW]=1;
	      // When elastic modulues increases, make desired area slightly smaller so cell does not bounce back too quick
	       c_A0[i] = 0.5 * cell_volume;
	       c_A0[cNEW] = 0.5 * cell_volume;
               c_growth[i]=0; c_growth[cNEW]=0;
	    }
        }
    }
}
//****************************************************************************
