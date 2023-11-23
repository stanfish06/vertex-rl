//****************************************************************************
//****************************************************************************
//***********************INITIAL STRUCTURE************************************
//****************************************************************************
//****************************************************************************
#include <iostream>
#include <fstream>
#include <sstream>

void set_initial_fromFile(){
    
    int *sides; sides = new int[13];
    int nrSides, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12;
    
    //FILE
    char filename2[100];
    snprintf(filename2, sizeof(char) * 100, "./initial/new_324_3.66_9.vt2d");
    FILE *file1; file1 = fopen(filename2, "rt");
    
    
    //V,E,P,C
    int nrV=0, nrE=0, nrC=0;
    fscanf(file1, "%d  %d  %d\n", &nrV, &nrE, &nrC);
    
    //perioXYZ
    fscanf(file1, "%lf  %lf\n", &perioXYZ[1], &perioXYZ[2]);
    
    //VERTICES
    float xx, yy, blank1, blank2;
    for(int i=1; i<=nrV; i++){
        fscanf(file1, "%lf  %lf  %lf  %lf\n", &xx, &yy, &blank1, &blank2);
        make_vertex(xx,yy);
    }

    //EDGES
    int v1, v2;
    float blank3;
    for(int i=1; i<=nrE; i++){
        fscanf(file1, "%d  %d  %lf\n", &v1, &v2, &blank3);
        make_edge(v1,v2,0);
    }

    //POLYGONS
    for(int i=0; i<=12; i++) sides[i]=0;
    for(int i=1; i<=nrC; i++){
        fscanf(file1, "%d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d\n", &nrSides, &b1, &b2, &b3, &b4, &b5, &b6, &b7, &b8, &b9, &b10, &b11, &b12);
        if(b1!=0) sides[1]=b1;
        if(b2!=0) sides[2]=b2;
        if(b3!=0) sides[3]=b3;
        if(b4!=0) sides[4]=b4;
        if(b5!=0) sides[5]=b5;
        if(b6!=0) sides[6]=b6;
        if(b7!=0) sides[7]=b7;
        if(b8!=0) sides[8]=b8;
        if(b9!=0) sides[9]=b9;
        if(b10!=0) sides[10]=b10;
        if(b11!=0) sides[11]=b11;
        if(b12!=0) sides[12]=b12;
        cell_sides(nrSides,sides,0);
    }
    delete [] sides;
    
    //MAKE CELLS
    for(int i=1; i<=Nc; i++) if(c[i][0]!=0) make_cell(i);
    
    //printf("%g  %g\n", perioXYZ[1], perioXYZ[2]);
    
    fclose(file1);
}
//****************************************************************************
int read_and_set_domain() {
    std::ifstream domain_coords("./initial/domain_coords.txt");
    std::ifstream domain_v_ns("./initial/domain_v_ns.txt");
    std::ifstream domain_v_ts("./initial/domain_v_ts.txt");

    if (!domain_coords.is_open()) {
        std::cerr << "Error opening the file" << std::endl;
	return 1;
    }

    int numLines = 0;
    std::string line;
    while (std::getline(domain_coords, line)) {
    	numLines++;
    }

    domain_coords.close();
    domain_coords.open("./initial/domain_coords.txt");

    v_D = new float*[numLines + 1];
    e_D = new int*[numLines + 1];
    e_normal_unit_vector = new float*[numLines + 1];
    e_tangent_unit_vector = new float*[numLines + 1];
    v_D[0] = new float[3];
    v_D[0][0] = numLines;
    e_D[0] = new int[2];
    e_D[0][0] = numLines;
    e_normal_unit_vector[0] = new float[3];
    e_normal_unit_vector[0][0] = numLines;
    e_tangent_unit_vector[0] = new float[3];
    e_tangent_unit_vector[0][0] = numLines;

    for (int i = 1; i <= numLines; i++) {
    	v_D[i] = new float[3];
	if (i == numLines) {
	    e_D[i] = new int[2];
	    e_D[i][0] = numLines;
	    e_D[i][1] = 1; 
	} else {
	    e_D[i] = new int[2];
	    e_D[i][0] = i;
	    e_D[i][1] = i + 1;
	}
	e_normal_unit_vector[i] = new float[3];
	e_tangent_unit_vector[i] = new float[3];
    }

    for (int i = 1; i <= numLines; i++) {
        std::getline(domain_coords, line);
	std::istringstream iss_v_D(line);
        float x, y;
	char comma;
        iss_v_D >> x >> comma >> y;
	//printf("%f, %f\n", x, y);
	v_D[i][1] = x;
	v_D[i][2] = y;

	std::getline(domain_v_ns, line);
        std::istringstream iss_e_n(line);
        iss_e_n >> x >> comma >> y;
	//printf("%f, %f\n", x, y);
        e_normal_unit_vector[i][1] = x;
        e_normal_unit_vector[i][2] = y;

	std::getline(domain_v_ts, line);
        std::istringstream iss_e_t(line);
        iss_e_t >> x >> comma >> y;
	//printf("%f, %f\n", x, y);
        e_tangent_unit_vector[i][1] = x;
        e_tangent_unit_vector[i][2] = y;	
    }

    return 0;
}

int read_and_set_interface() {
    std::ifstream interface_coords("./initial/interface_coords.txt");
    std::ifstream interface_v_ns("./initial/interface_v_ns.txt");
    
    if (!interface_coords.is_open()) {
        std::cerr << "Error opening the file" << std::endl;
        return 1;
    }
    int numLines = 0;
    std::string line;
    while (std::getline(interface_v_ns, line)) {
        numLines++;
    }

    interface_v_ns.close();
    interface_v_ns.open("./initial/interface_v_ns.txt");
    
    v_interface = new float*[numLines + 1];
    e_normal_unit_vector_interface = new float*[numLines + 1];
    v_interface[0] = new float[3];
    v_interface[0][0] = numLines;
    e_normal_unit_vector_interface[0] = new float[3];
    e_normal_unit_vector_interface[0][0] = numLines;
 
    for (int i = 1; i <= numLines; i++) {
    	v_interface[i] = new float[3];
	e_normal_unit_vector_interface[i] = new float[3];
    }

    std::string line_coords;
    std::string line_v_ns;
    std::getline(interface_coords, line_coords);
    std::istringstream iss_v_interface(line_coords); 
    for (int i = 1; i <= numLines; i++) {
        float x, y;
	char comma;
	if (i == 1) {
            iss_v_interface >> x >> comma >> y;
	} else {
	    iss_v_interface >> comma >> x >> comma >> y;
	}
	//printf("%f, %f\n", x, y);
	v_interface[i][1] = x;
	v_interface[i][2] = y;

	std::getline(interface_v_ns, line_v_ns);
        std::istringstream iss_e_n(line_v_ns);
        iss_e_n >> x >> comma >> y;
	//printf("%f, %f\n", x, y);
        e_normal_unit_vector_interface[i][1] = x;
        e_normal_unit_vector_interface[i][2] = y;
    }
    
    return 0;
}

void set_initial_singleCell(){
    
    int *sides; sides = new int[101];
    
    float min_x = 0;
    float max_x = -1;
    float min_y = 0;
    float max_y = -1;
    
    for (int i = 1; i <= v_D[0][0]; i++) {
        float x = v_D[i][1];
	float y = v_D[i][2];

	if (x > max_x) {
	    max_x = x;
	}

	if (y > max_y) {
	    max_y = y;
	}
    } 
    // sampling box
    //float width = scale_factor * (4*sqrt(2) + 4);
    //float height = scale_factor * (14 + 4*sqrt(2));
    float width = max_x;
    float height = max_y; 
    // number of random samples
    srand(646);
    int nrC_try =300;
    int nrC_per_row = int(width);
    int nrC_per_col = int(height);
    //V,E,P,C
    int nrV=max_num_edges, nrE=max_num_edges, nrC=0;
    //int nrV=12, nrE=12, nrC=1;

    //perioXYZ
    perioXYZ[1]=L;
    perioXYZ[2]=L;
    
    //VERTICES
    float R0=0.3;
    //for(int i=1; i<=nrV; i++) make_vertex(R0*cos(2*M_PI*i/nrV+0.1)+0.5*perioXYZ[1],R0*sin(2*M_PI*i/nrV+0.1)+0.5*perioXYZ[2]-L/3.);
   // compute_normal_vector();
    print_normal_vector();
   // int c_y = nrC_per_col * 2 - 1;
    /*
    float grid_width = width / nrC_per_row;
    float grid_height = height / nrC_per_col;
    for (int i = 0; i < nrC_per_row; i++) {
        for (int j = 0; j < nrC_per_col; j++) {
	    float x = i * grid_width;
	    float y = j * grid_height;
	    float c_x = x + rnd() * grid_width;
	    float c_y = y + rnd() * grid_height;
	    if (in_domain_some_distance_to_edge(c_x, c_y)) {
                for(int i=1; i<=nrV; i++) make_vertex(R0*cos(2*M_PI*i/nrV+0.1)+c_x,R0*sin(2*M_PI*i/nrV+0.1)+c_y);
                nrC++;
            }
	}
    }
    */
    
    for(int n = 1; n <= nrC_try; n++) {
        float c_x = rnd() * width;
        float c_y = rnd() * height;
	//c_x = 32.860644;
	//c_y =	5.582756;
//	c_y = 3.186068;
	//int c_x = n % nrC_per_row * 2;
	//if (c_x == 0) {
	//   c_y = c_y - 2;
	//} 
	//float c_x_jit = c_x*1.0 + sin(2.0 * M_PI * rnd());
        //float c_y_jit = c_y*1.0 + sin(2.0 * M_PI * rnd()); 
	// make sure to not random dropping in proliferative zone
        if (in_domain_some_distance_to_edge(c_x, c_y)) {
	    //printf(" Yes\n");
	    //std::cout << in_domain_some_distance_to_edge(c_x, c_y) << std::endl;
            for(int i=1; i<=nrV; i++) make_vertex(R0*cos(2*M_PI*i/nrV+0.1)+c_x,R0*sin(2*M_PI*i/nrV+0.1)+c_y);
            nrC++;
        }
    }
    
    

    //for(int i=1; i<=nrV; i++) make_vertex(R0*cos(2*M_PI*i/nrV+0.1)+ 2*sqrt(2) + 0.5, R0*sin(2*M_PI*i/nrV+0.1)+ 10);
    //for(int i=1; i<=nrV; i++) make_vertex(R0*cos(2*M_PI*i/nrV+0.1)+ 2*sqrt(2) + 2, R0*sin(2*M_PI*i/nrV+0.1)+ 10);
    //for(int i=1; i<=nrV; i++) make_vertex(R0*cos(2*M_PI*i/nrV+0.1)+ 2*sqrt(2) + 3.5, R0*sin(2*M_PI*i/nrV+0.1)+ 10);
    //for(int i=1; i<=nrV; i++) make_vertex(R0*cos(2*M_PI*i/nrV+0.1)+5,R0*sin(2*M_PI*i/nrV+0.1)+7);
    //for(int i=1; i<=nrV; i++) make_vertex(R0*cos(2*M_PI*i/nrV+0.1)+3,R0*sin(2*M_PI*i/nrV+0.1)+7);
    //for(int i=1; i<=nrV; i++) make_vertex(R0*cos(2*M_PI*i/nrV+0.1)+4,R0*sin(2*M_PI*i/nrV+0.1)+7);
/*
    make_vertex(13.2093, 21.5031);
    make_vertex(14.0603, 22.2015);
    make_vertex(13.1219, 23.5359);
    make_vertex(12.3060, 23.3271);
    make_vertex(12.0049, 23.0383);
    make_vertex(11.8107, 22.5001);
    make_vertex(11.7645, 22.4839);
    make_vertex(10.4778, 21.7091);
    make_vertex(10.3984, 22.6628);
    make_vertex(11.0443, 23.3022);
    make_vertex(11.96, 23.2568);
    make_vertex(11.8835, 22.4964);
*/	
    //EDGES
    int v1, v2;
    for(int k = 0; k < nrC; k++) {
        for(int i=1 + nrV * k; i<nrE  + nrV * k; i++) make_edge(i,i+1,0);
        make_edge(nrE + nrV * k, 1 + nrV * k, 0);
    }

    //POLYGONS
    for(int m=0; m<nrC; m++) {
        for(int i=0; i<=100; i++) sides[i]=0;
        for(int i=1; i<=nrE; i++) sides[i]=i + m * nrE;
        cell_sides(nrE,sides,0);
    }
    //for(int i=0; i<=100; i++) sides[i]=0;
    //for(int i=1; i<=nrE; i++) sides[i]=i;
    //cell_sides(nrE,sides,0);
    delete [] sides;
    
    //MAKE CELLS
    printf("%d", Nc);

    for(int i=1; i<=Nc; i++) if(c[i][0]!=0) {
        make_cell(i);
        calc_cell_center(i);
    }

    //Domain polygon
    //v_D = {{6, 0, 0}, {0, 1, 2}, {0, 4, 1}, {0, 8, 2}, {0, 8, 6}, {0, 4, 12}, {0, 1, 8}};
    //e_D = {{6, 0}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 1}};
    //e_normal_unit_vector = {{6, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    //printf(e_normal_unit_vector);
    //UPDATE
    calc_edges();
    update_shoulder_vertices();
}
//****************************************************************************
void set_initial_triangles(int NrTri, float lsmall){
    
    float **centers;
    centers = new float*[array_max+1];
    for(int i=0; i<=array_max; i++) centers[i] = new float[3];
    for(int i=0; i<=array_max; i++) for(int j=0; j<=2; j++) centers[i][j]=0;
    
    //V,E,P,C
    int nrV=3*NrTri, nrE=3*NrTri, nrC=NrTri;
    
    //perioXYZ
    perioXYZ[1]=sqrt(1.*NrTri);
    perioXYZ[2]=sqrt(1.*NrTri);
    
    for(int i=1; i<=nrC; i++){
        
        float xc, yc;
        while(1<2){
            xc=perioXYZ[1]*rnd(); yc=perioXYZ[2]*rnd();
            int check=-1;
            for(int j=1; j<i; j++){
                float xr=centers[j][1]; float yr=centers[j][2];
                float dd=torus_distance(xc,yc,xr,yr);
                if(dd<0.2) { check=1; break; }
            }
            if(check==-1) break;
        }
        centers[i][2]=yc;//check overlaps
        centers[i][1]=xc;//check overlaps
        
        //vertices
        int v1=make_vertex(xc-lsmall,yc-lsmall); torus_vertex(v1);
        int v2=make_vertex(xc+lsmall,yc-lsmall); torus_vertex(v2);
        int v3=make_vertex(xc,yc+2*lsmall); torus_vertex(v3);
        
        //edges
        int e1=make_edge(v1,v2,0);
        int e2=make_edge(v2,v3,0);
        int e3=make_edge(v3,v1,0);
        
        //polygon
        int *sides; sides = new int[101]; for(int i=0; i<=100; i++) sides[i]=0;
        sides[1]=e1;
        sides[2]=e2;
        sides[3]=e3;
        cell_sides(3,sides,0);
        delete [] sides;
    }
    
    //MAKE CELLS
    for(int i=1; i<=Nc; i++) if(c[i][0]!=0) make_cell(i);
    
    printf("%g  %g\n", perioXYZ[1], perioXYZ[2]);
    
    for(int i=0; i<=array_max; i++) delete [] centers[i];
    delete [] centers;
}
//****************************************************************************
