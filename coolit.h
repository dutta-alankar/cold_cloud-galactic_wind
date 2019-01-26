#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MWIDTH  2
#define MHEIGHT 1001

double lambda(double T){

    FILE* fp;
    int height, width, ii, jj;
    float array[MHEIGHT][MWIDTH];

    if((fp = fopen("cooltable.dat", "r")) == NULL)
        exit(1);

    for(jj=0; jj<MHEIGHT; jj++)
        for(ii=0; ii<MWIDTH; ii++)
            if(fscanf(fp, "%e", &array[jj][ii])!=1) exit(1);
                
    fclose(fp);
/*
    for(jj=0; jj<MHEIGHT; jj++){
        for(ii=0; ii<MWIDTH; ii++)
            printf ("%e    ", array[jj][ii]);
        printf("\n");
    }
*/



    for(jj=1; jj<MHEIGHT; jj++){
        if ((T-array[jj][0])<0 ){ 
          if (abs(T-array[jj-1][0])>(T-array[jj][0])) return array[jj][1];
          else return array[jj-1][1];
        }
    }
    return array[MHEIGHT-1][1];
}
