//***************************************************************
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
//***************************************************************
int main (int argc, char *argv[])
{

    double p0=3.82;
    char fileName[50];
    
    //MASTER
    FILE *outfile;
    snprintf(fileName, sizeof(char) * 50, "master_%g.slurm", p0);
    outfile = fopen(fileName, "at");
    for(int par1=0; par1<=50; par1++) for(int par2=0; par2<=50; par2++) fprintf(outfile, "sbatch r_%g_%d_%d.slurm\n", p0, par1, par2);
    fclose(outfile);
    
    
    
    //INDIVIDUAL
    for(int par1=0; par1<=50; par1++){
        for(int par2=0; par2<=50; par2++){
        	snprintf(fileName, sizeof(char) * 50, "r_%g_%d_%d.slurm", p0, par1, par2);
        	outfile = fopen(fileName, "wt");
        	fprintf(outfile, "#!/bin/bash\n\n");
        	fprintf(outfile, "#SBATCH --nodes=1\n");
        	fprintf(outfile, "#SBATCH --ntasks-per-node=1\n");
        	fprintf(outfile, "#SBATCH --mem=4000\n");
        	fprintf(outfile, "#SBATCH --time=05:00:00\n\n");
        	fprintf(outfile, "time ./main.o %g %d %d\n", p0, par1, par2);
        	fclose(outfile);
    	}
    }
    
    return 0;
}
//***************************************************************
