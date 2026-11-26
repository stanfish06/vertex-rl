#include <stdio.h>
#include <cstring>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <time.h>

//VERTICES
float **v;
float **vS;
//.................
int **v_edges;
int **v_cells;
//.................
float **v_F;

//EDGES
int **e;
//.................
int **e_cells;
int **e_Vcells;
//.................
float *e_g;
float *e_length;
float *e_lengthS;
float *e_dl;
float *e_clock;
float **e_t;
float **e_n;
float **e_center;

//CELLS
int **c;
int *c_dead;
int *c_shrink;
//.................
int **c_vertices;
int **c_Vedges;
//.................
float *c_area;
float *c_A0;
float **c_center;
int *c_growth;

//PURSE STRING
float eFg=0;
float eGamma=0;

//MISC
bool force_topology_changed = true;
int **edge_vertex_int; // indexed by edge, then vertex
float **ve_next_time;
float **vv_next_time;
int Nv, Ne, Nc;
float *perioXYZ;
#ifndef VERTEX_ARRAY_CAPACITY
#define VERTEX_ARRAY_CAPACITY 50000
#endif
const int array_max = VERTEX_ARRAY_CAPACITY;
static_assert(array_max > 0, "vertex array capacity must be positive");
int **ids, **exist;
float wA, wL, wS;

//Domain
//float scale_factor = 4.0;
/*
float v_D[11][3] = {{10, 0, 0}, {0, scale_factor * 2*sqrt(2), 0}, {0, scale_factor * 2*sqrt(2), scale_factor * 10}, {0, 0, scale_factor * (10 + 2*sqrt(2))}, {0, 0, scale_factor * (14 + 2*sqrt(2))}, {0, scale_factor * 2*sqrt(2), scale_factor * (14 + 2*2*sqrt(2))}, {0, scale_factor * (2*sqrt(2) + 4), scale_factor * (14 + 2*2*sqrt(2))}, {0, scale_factor * (4*sqrt(2) + 4), scale_factor * (14 + 2*sqrt(2))}, {0, scale_factor * (4*sqrt(2) + 4), scale_factor * (10 + 2*sqrt(2))}, {0, scale_factor * (2*sqrt(2) + 4), scale_factor * 10}, {0, scale_factor * (2*sqrt(2) + 4), 0}};
int e_D[11][2] = {{10, 0}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 9}, {9, 10}, {10, 1}};

int e_D_head[8] = {7, 2, 3, 4, 5, 6, 7, 8};
int e_D_tail[4] = {3, 1, 9, 10};

// normal vectors point into the domain
float e_normal_unit_vector[11][3] = {{10, 0, 0}, {0, 1, 0}, {0, 0.707107, 0.707107}, {0, 1, 0}, {0, 0.707107, -0.707107}, {0, 0, -1}, {0, -0.707107, -0.707107}, {0, -1, 0}, {0, -0.707107, 0.707107}, {0, -1, 0}, {0, 0, 1}};
*/
float **v_D;
int **e_D;
float **e_normal_unit_vector;
float **e_tangent_unit_vector;

float **e_normal_unit_vector_interface;
float **v_interface;
// Angles and direction of division
int lateral_normal[2] = {1, 9};
float boundary_x[2] = {(2*sqrt(2)) * scale_factor, (2*sqrt(2) + 4) * scale_factor};

#include "_torus.h"
#include "_arrays.h"
#include "_free_id.h"
#include "_dissolve.h"
#include "_domain_initialization.h"
#include "_domain_constraint.h"
#include "_allocate.h"
#include "_subdomain_determination.h"
#include "_create.h"
#include "_rndom.h"
#include "_edge_division.h"
#include "_edge_merge.h"
#include "_crossing_edge_resolution.h"
#include "_edge_manipulation.h"
#include "_angle_of_division.h"
#include "_cell_division.h"
#include "_edge.h"
#include "_output.h"
#include "_force_workspace.h"
#include "_force_length.h"
#include "_force_area.h"
#include "_force_steric.h"
#include "_force_shoulder.h"
#include "_equation_of_motion.h"
#include "_initial_structure.h"
