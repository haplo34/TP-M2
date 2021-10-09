#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        printf("%s: oops - no command line argument given\n", argv[0]);
        return -1; /* error */
    }

    double a = strtod(argv[1], NULL);
    double b = strtod(argv[2], NULL);
    double N = strtod(argv[3], NULL);
    
    for(double i=0; i<N; i++)
    {
        double c = a/b;
    }
    
    printf("arguments supplied are: %s, %s, %s\n", argv[1], argv[2], argv[3]);
    return 0;
   
}
