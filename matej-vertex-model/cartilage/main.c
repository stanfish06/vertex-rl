#include <cmath>
#include <stdio.h>
#include <readline/readline.h>
#include <readline/history.h>

//****************INPUT PARAMETERS************************
//MODEL PARAMETERS
float k0=99;      //cell-cell repulsion
float kS=2;      //ECM shoulder repulsion
float lS=0.5;      //ECM shoulder thickness
float kC=0.5;
//float e_g0=0.25;    //cell line tension (should be below 0.434313)
float Atarget=1;   //cell target area (should be kept 1=normalized)
float kA=1.2; //1.2       //cell-area elastic modulus (should be kept 1=normalized)
float e_g0=0.7; //0.4    //cell line tension (should be below 0.434313)
// Threshold of line tension should be calculated by balancing with area compression, slightly higher than threshold is fine as desired area will change.
// However, if initial area is too small and cell's desired area changes slow, cell area might shrink to 0 as there is no nonzero energy minimum 
float tMAX=2000;   //total integration time
float outFreq=0.5;   //frequency of outputing
float h0=0.001;    //time step (should not be too big)
float L=4;        //width of the simulation box

//PARAMETERS FOR EDGE MANIPULATION
float lthUP=0.2;   //if edge longer than this, it gets divided
float lthDOWN=0.01;//if edge shorter than this, its vertices are merged
float lALERT=0.4;  //width of the region around edges, where vertex-edge interactions are computed
float lDEEP=0.5;
//OTHER PARAMETERS
float h;
float Time;
float maxDR=0.001;
float ddTH, ddTHS, vMAX;
float TimeCount=0;
int kount=0;

// Parallelism parameters
int num_threads_count = 16;
int num_threads_division = 16;
int num_threads_edge_manipulation = 10;
int num_threads_force = 16;
int num_threads_calculation = 16;
int num_threads_update = 16;

float scale_factor = 2.0;

// Zone interfaces
float proliferative_zone_interface[] = {7 * scale_factor, 11 * scale_factor};
float hypertrophic_zone_interface[] = {4 * scale_factor, 7 * scale_factor};
float resting_zone_interface[] = {11 * scale_factor, (14 + 2 * 2 * sqrt(2)) * scale_factor};
float ossification_interface[] = {0 * scale_factor, 4 * scale_factor};
float pericondrium_thickness = 1.0;

float proliferative_zone_area = (proliferative_zone_interface[1] - proliferative_zone_interface[0]) * 4 * scale_factor;
// Polygon interface (all subpolygons have to be convex)
float polygon_interface_y = 10 * scale_factor;

// Region ECM (scale factor for lS)
float proliferative_ECM = 2;
float hypertrophic_ECM = 3;
float pericondrium_ECM = 1;
float ossification_ECM = 3;
float resting_ECM = 1;

float ECM[] = {
    resting_ECM,
    proliferative_ECM,
    hypertrophic_ECM,
    ossification_ECM,
    pericondrium_ECM
};

// Core (scale factor for lALERT)
float proliferative_CORE = sqrt(1.5);
float hypertrophic_CORE = sqrt(3);
float pericondrium_CORE = 1;
float ossification_CORE = sqrt(3);
float resting_CORE = 1;

float CORE[] = {
    resting_CORE,
    proliferative_CORE,
    hypertrophic_CORE,
    ossification_CORE,
    pericondrium_CORE
};


// Region cell_volume (scale factor for lS)
float proliferative_CV = 1.5;
float hypertrophic_CV = 2.5;
float pericondrium_CV = 1;
float ossification_CV = 2.5;
float resting_CV = 1;

float CV[] = {
    resting_CV,
    proliferative_CV,
    hypertrophic_CV,
    ossification_CV,
    pericondrium_CV
};

// Region line tension (scale factor for eg0)
float proliferative_EG = 0.75;
float hypertrophic_EG = 1;
float pericondrium_EG = 1;
float ossification_EG = 1;
float resting_EG = 1;

float EG[] = {
    resting_EG,
    proliferative_EG,
    hypertrophic_EG,
    ossification_EG,
    pericondrium_EG
};



// Edge threshold
int max_num_edges = 14;
int min_division_n_edges = 12;
int min_exist_n_edges = 6;
//*******************DECLARATIONS*************************
#include "_functions.h"
#include <omp.h>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
//*******************FUNCTION****************************
void parseAndSet(const std::string &input) {
    std::istringstream iss(input);
    std::string command, parameter;
    float value;
    iss >> command;
    if (command == "set") {
        iss >> parameter >> value;
	//std::cout << parameter << std::endl;
	//printf("value %f\n", value);
	//std::cout << std::to_string(parameter == "k0") << std::endl;
        if (parameter  == "k0") k0 = value;
	else if (parameter  == "kS") kS = value;
	else if (parameter  == "lS") lS = value; 
	else if (parameter  == "kA") kA = value;	
	else if (parameter  == "e_g0") e_g0 = value;
	else if (parameter  == "tMAX") tMAX = value;
	else if (parameter  == "outFreq") outFreq = value;
	else if (parameter  == "maxDR") maxDR = value;
	else if (parameter  == "lthUP") lthUP = value;
	else if (parameter  == "lthDOWN") lthDOWN = value;
	else if (parameter  == "lALERT") lALERT = value;
	else if (parameter  == "max_num_edges") max_num_edges = value;
	else if (parameter  == "min_division_n_edges") min_division_n_edges = value;
	else if (parameter  == "min_exist_n_edges") min_exist_n_edges = value;
	else  std::cout << "Unknown parameter: " << parameter << std::endl;	
	std::cout << std::endl;
	return;
    }	
    return;
}

void parameter_list() {        
printf(
    " Name                           Description                     Current value   \n"
    "  1. k0                         Description of parameter 1      %f              \n"
    "  2. kS                         Description of parameter 2      %f              \n"
    "  3. lS                         Description of parameter 3      %f              \n"
    "  4. kA                         Description of parameter 4      %f              \n"
    "  5. e_g0                       Description of parameter 5      %f              \n"
    "  6. tMAX                       Description of parameter 6      %f              \n"
    "  7. outFreq                    D                               %f              \n"
    "  8. maxDR                      D                               %f              \n"
    "  9. lthUP                      D                               %f              \n"
    " 10. lthDOWN                    D                               %f              \n"
    " 11. lALERT                     D                               %f              \n"
    " 12. max_num_edges              D                               %d              \n"
    " 13. min_division_n_edges       D                               %d              \n"
    " 14. min_exist_n_edges          D                               %d              \n",
    k0, kS, lS, kA, e_g0, tMAX, outFreq, maxDR, lthUP, lthDOWN, lALERT, max_num_edges, min_division_n_edges, min_exist_n_edges
);
}

void domain_setting() {
    if (v_D == nullptr) {
    	printf("Domain not set yet. Run specify domain first.\n");
	return;
    }
    std::cout << "Domain: " << std::endl;
    for (int i = 1; i <= v_D[0][0]; i++) {
	std::cout << "Domain edge " << i << " :" << std::endl;
	printf("Edge coords: (%f, %f)--(%f, %f)\n", v_D[e_D[i][0]][1], v_D[e_D[i][0]][2], v_D[e_D[i][1]][1], v_D[e_D[i][1]][2]);
	printf("Edge normal vectors: (%f, %f)\n", e_normal_unit_vector[e_D[i][0]][1], e_normal_unit_vector[e_D[i][0]][2]);
	printf("Edge tangent vectors: (%f, %f)\n", e_tangent_unit_vector[e_D[i][0]][1], e_tangent_unit_vector[e_D[i][0]][2]);
        std::cout << std::endl;
    }
    std::cout << std::endl;    
}

void interface_setting() {
    if (v_interface == nullptr) {
        printf("Interface not set yet. Run specify interface first.\n");
        return;
    }
    std::cout << "Interface: " << std::endl;
    for (int i = 1; i <= v_interface[0][0]; i++) {
        std::cout << "Interface " << i << " :" << std::endl;
        printf("Interface coords: (%f, %f)\n", v_interface[i][1], v_interface[i][2]);
        printf("Interface normal vectors: (%f, %f)\n", e_normal_unit_vector_interface[i][1], e_normal_unit_vector_interface[i][2]);
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

//********************MAIN********************************
int main(int argc, char *argv[]){
    const char* helpManual = R"(
Usage: 2D vertex model for cartilage simulation

Commands:
  help               			Display this help message
  ls parameters      			List all parameters
  set [s]            			Set value for parameter s. Use <ls parameters> to see full list of parameters.
  specify domain     			Specify the domain of simulation. Remember to enter points in clockwise order. Don't enter same point twice. 
  ls domain				Print current domain setting.
  specify interface                     Specify the interface of simulation. Remember to enter points in clockwise order. Don't enter same point twice.
  ls interface                          Print current interface setting.
  run		     			Run simulation with current setting.    

Examples:
)";
const char* parameterManual = R"(
  Name                          Description			Current value
  1. k0          		Description of parameter 1	
  2. kS         		Description of parameter 2
  3. lS  			Description of parameter 3
  4. kA  			Description of parameter 4 
  5. e_g0          		Description of parameter 5
  6. tMAX                       Description of parameter 5
  7. outFreq                    D
  8. maxDR			D
  9. lthUP			D
 10. lthDOWN			D
 11. lALERT			D
 12. max_num_edges		D
 13. min_division_n_edges	D
 14. min_exist_n_edges		D
)";
    printf("Press Ctrl+c or enter exit to quit\n");
    char* buf;
    while ((buf = readline(">> ")) != nullptr) {
	std::string buf_str = buf;
        if (strcmp(buf, "help") == 0) {
	    std::cout << helpManual << std::endl;
	} else if (strcmp(buf, "exit") == 0) {
	    exit(0);
	} else if (strcmp(buf, "ls parameters") == 0) {
 	    //std::cout << parameterManual << std::endl;
	    parameter_list();
    	} else if (buf_str.substr(0, 3) == "set") {
	    parseAndSet(buf_str);
	} else if (strcmp(buf, "specify domain") == 0) {
            if (read_and_set_domain() == 1) {
	    	printf("Domain setting fails.\n");
	    }
	} else if (strcmp(buf, "ls domain") == 0) {
	    domain_setting();
	} else if (strcmp(buf, "specify interface") == 0) {
	    if (read_and_set_interface() == 1) {
                printf("Interface setting fails.\n");
            }
        } else if (strcmp(buf, "ls interface") == 0) {
            interface_setting();
        } else if (strcmp(buf, "run") == 0) {
            printf("Starting...\n");
	    goto start_simulation;
	} else {
	    printf("Unknown command\n");
	}
//	free(buf);
    }   
    return 0;
   
start_simulation:    
    //INITIALIZE
    ddTH=sqrt(pow(lthUP/2,2)+pow(lALERT,2));
    ddTHS=sqrt(pow(lthUP/2,2)+pow(lS+lALERT,2));
    vMAX=2.*maxDR/h0;
    allocate();
    if (read_and_set_domain() == 1) {printf("Domain setting fails"); return 1;}
    if (read_and_set_interface() == 1) {printf("Interface setting fails"); return 1;}
   // double test_x = 37.710607;
   // double test_y = 5.487489;
   // std::cout << in_domain_some_distance_to_edge(37.710607, 5.487489) << std::endl;
   // std::cout << in_domain_some_distance_to_edge(test_x, test_y) << std::endl;
    // return 1;
    set_initial_singleCell();
    //return 1;
    //*************************
    int MAX = 20; 
    int count = 0;
    //#pragma omp parallel num_threads(MAX)
    //{
    //   #pragma omp atomic
    //   count++;
    //}
    //printf("Number of threads: %d\n", count); 
    //TIME LOOP
    Time=0;
   // kount++;
    //out_tissue(kount);
   // return 1;
    while(Time<tMAX){
	//printf("Number of threads: %d\n", count);
        //GROWTH & DIVISION
        //for (int i = 1; i <= Nc; i++) {
        //    resolve_cross_cell(i);
        //}

	//grows and divides cells according to custom rules
        growth_division();
        
        //EQUATION OF MOTION
        //calculates forces on vertices and performs 1 integration step
        if(eqOfMotion()==0) break;
        
        //TOPOLOGICAL TRANSITIONS
        //makes sure that edges are not too long or too short
        
	//edge_manipulation();
        //printf("%d", Nc);
	//error_vertex_deletion();
	for (int i = 1; i <= Nc; i++) {
	    resolve_cross_cell(i);
	}
        edge_manipulation();
    }
  
    //DEALLOCATE
    deallocate();
    //*************************

    return 0;
}
//********************************************************
