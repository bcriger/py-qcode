#include "squoct_dist.h"
#include "stdio.h"

int main(int argc, char const *argv[])
{
    unsigned short int x1, y1, x2, y2, sz_x, sz_y;
    num_t d;
    char synd_type, quit_query;
    while(1){
        
        puts("Enter x1"); scanf("%hu", &x1);
        puts("Enter y1"); scanf("%hu", &y1);
        puts("Enter x2"); scanf("%hu", &x2);
        puts("Enter y2"); scanf("%hu", &y2);
        puts("Enter sz_x"); scanf("%hu", &sz_x);
        puts("Enter sz_y"); scanf("%hu", &sz_y);
        puts("Enter synd_type"); scanf(" %c", &synd_type);

        printf("First Point: (%d, %d)\n", x1, y1);
        printf("Second Point: (%d, %d)\n", x2, y2);
        printf("Syndrome Type: %c\n", synd_type);

        d = dist((num_t)x1, (num_t)y1, (num_t)x2, (num_t)y2, 
                    (num_t)sz_x, (num_t)sz_y, synd_type);
        
        printf("Resulting Distance: %hu\n", (unsigned short int)d);
        
        puts("Enter q to quit");
        scanf(" %c", &quit_query);
        if (quit_query=='q') break;
    }
    return 0;
}