//****************************************************************************
//****************************************************************************
//****************************T1 TRANSITION***********************************
//****************************************************************************
//****************************************************************************
//#include <omp.h>
int edge_manipulation()
{
//#pragma omp parallel for num_threads(num_threads_edge_manipulation)
    for(int i=1; i<=Ne; i++) if(exist[2][i]==1 && e_clock[i]==1){
        if(e_length[i]>lthUP && c[e_cells[i][1]][0] < max_num_edges) EDGE_DIVISION(i);
        else if(e_length[i]<lthDOWN && c[e_cells[i][1]][0]>min_exist_n_edges) {EDGE_MERGE(i);} //printf("oh no!\n"); }
    }
    for(int i=1; i<=Ne; i++) if(exist[2][i]==1) e_clock[i]=1;
    return 0;
}
//****************************************************************************
