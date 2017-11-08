#include <stdio.h>
#include <stdlib.h>

#include "besselI.h"

int main(int argc, char *argv[])
{
    if(argc < 3)
    {
        fprintf(stderr, "%s order argument\n", argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    double x = atof(argv[2]);

    printf("%.15g\n", besselI(n,x));
}
