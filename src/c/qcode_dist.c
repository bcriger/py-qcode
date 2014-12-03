#include "qcode_dist.h"

#ifdef _DEBUG
#include "stdio.h"
#endif
#ifdef _DEBUG_NEIGHBOURS
#include "stdio.h"
#endif

num_t i_min_2(num_t a, num_t b){
    if (a < b) return a;
    else return b;
}

num_t i_min_4(num_t a, num_t b, num_t c, num_t d){
    num_t ret_val = a;
    if (b < ret_val) ret_val = b;
    if (c < ret_val) ret_val = c;
    if (d < ret_val) ret_val = d;
    return ret_val;
}

num_t min_abs_sub(num_t x1, num_t x2, num_t sz_x){
    /*
    Special care must be taken here to avoid unsigned integer wrapping.
    */
    int abs_sub = 0;
    
    if (x1 > x2) abs_sub = x1 - x2;
    else abs_sub = x2 - x1;
    
    #ifdef _DEBUG
    puts("In min_abs_sub, evaluating absolute difference of two 1D coordinates:");
    printf("x1, x2, diff: %hu, %hu, %hu \n", 
        (unsigned short int) x1,
        (unsigned short int) x2,
        (unsigned short int) abs_sub);
    puts("-------------------");
    #endif
    
    return i_min_2(abs_sub, sz_x - abs_sub);
}

bool is_sq_cent(num_t x, num_t y){
    /*
    Acts on coordinates in the virtual (square) lattice, returning true
    if the coordinates are 'skew' (not both even and not both odd).
    */
    if(x % 2 == y % 2) return false;
    else return true;
}
    
num_t octagonal_dist(num_t x1, num_t y1, num_t x2, num_t y2, num_t sz_x, num_t sz_y){
    /*
    The distance between two octagonal checks is given by the Manhattan
    norm after the co-ordinates have been shrunk down to the virtual 
    lattice. This function evaluates that distance.
    */
    int num_steps = 0;
    num_steps += min_abs_sub(x1, x2, sz_x);
    
    #ifdef _DEBUG
    puts("In octagonal_dist");
    printf("Distance along x-direction from %hu to %hu (wrapping at %hu): %hu\n",
            (unsigned short int)x1, (unsigned short int)x2,
            (unsigned short int)sz_x, (unsigned short int)num_steps);
    #endif
    
    num_steps += min_abs_sub(y1, y2, sz_y);
    
    #ifdef _DEBUG
    printf("Total distance adding y-data from %hu to %hu (wrapping at %hu): %hu\n",
            (unsigned short int)y1, (unsigned short int)y2,
            (unsigned short int)sz_y, (unsigned short int)num_steps);
    puts("-------------------");
    #endif

    return num_steps;
}

num_t squoct_dist(num_t x1, num_t y1, num_t x2, num_t y2, num_t sz_x, num_t sz_y, char synd_type){
    /*
    This function returns the distance on the dual lattice for the 
    concatenated toric/[[4,2,2]] code defined in 
    py_qcode/src/py_qcode/(lattice.py, code.py). It is fairly 
    complicated, so implementation in pure python uses more time
    than c++-implemented minimum-weight matching. This C 
    implementation will be compiled into a shared library and 
    accessed using ctypes in lattice.py.
    */

    //VARIABLE DECLARATIONS
    num_t num_steps = 0; //Return Value
    //positions neighbouring square coordinates
    num_t x1_u, x1_d, x2_u, x2_d, y1_u, y1_d, y2_u, y2_d;

    /*
    Step 0: Short circuit for points which are identical:
    */
    //
    if (x1 == x2 && y1 == y2)
    {
        return 0;
    }
    else
    {
    /*
    Step 1: Each of the coordinates of the dual lattice has been 
    subjected to an affine map x -> (3x + 1), so that:
    + the lowest represented coordinate value will be 0
    + no collisions will occur between the neighbourhoods of each point
      on this dual lattice.
    We invert this map:  
    */
    x1 = (x1 - 1) / 3; x2 = (x2 - 1) / 3; sz_x /= 3;
    y1 = (y1 - 1) / 3; y2 = (y2 - 1) / 3; sz_y /= 3;
    
    #ifdef _DEBUG
    puts("In dist, these coordinates are supposed to be on the virtual square lattice:");
    printf("x1, x2, sz_x: %hu, %hu, %hu\n", 
        (unsigned short int) x1, 
        (unsigned short int) x2,
        (unsigned short int) sz_x);
    printf("y1, y2, sz_y: %hu, %hu, %hu\n", 
        (unsigned short int) y1, 
        (unsigned short int) y2,
        (unsigned short int) sz_y);
    #endif

    /*
    Step two: We compute the neighbouring co-ordinates of all input coordinates,
    in case some of the inputs represent squares:
    */
    if (x1 % 2){ //x coord is odd
        if (synd_type == 'Z'){
            x1_u = (x1 + 1) % sz_x; x1_d = (x1 - 1);
            y1_u = y1; y1_d = y1; 
        }
        else if (synd_type == 'X'){
            x1_u = x1; x1_d = x1;
            y1_u = (y1 + 1) % sz_y; if (y1 == 0) y1_d = sz_y - 1; else y1_d = y1 - 1; 
        }
    }
    else { //x coord is even
        if (synd_type == 'Z'){
            y1_u = (y1 + 1) % sz_y; if (y1 == 0) y1_d = sz_y - 1; else y1_d = (y1 - 1);
            x1_u = x1; x1_d = x1; 
        }
        else if (synd_type == 'X'){
            y1_u = y1; y1_d = y1;
            x1_u = (x1 + 1) % sz_x; if (x1 == 0) x1_d = sz_x - 1; else x1_d = x1 - 1; 
        }
    }
    if (x2 % 2){ //x coord is odd
        if (synd_type == 'X'){
            y2_u = (y2 + 1) % sz_y; if (y2 == 0) y2_d = sz_y - 1; else y2_d = y2 - 1;
            x2_u = x2; x2_d = x2; 
        }
        else if (synd_type == 'Z'){
            y2_u = y2; y2_d = y2;
            x2_u = (x2 + 1) % sz_x; if (x2 == 0) x2_d = sz_x - 1; else x2_d = x2 - 1; 
        }
    }
    else { //x coord is even
        if (synd_type == 'Z'){
            y2_u = (y2 + 1) % sz_y; if (y2 == 0) y2_d = sz_y - 1; else y2_d = y2 - 1;
            x2_u = x2; x2_d = x2; 
        }
        else if (synd_type == 'X'){
            y2_u = y2; y2_d = y2;
            x2_u = (x2 + 1) % sz_x; if (x2 == 0) x2_d = sz_x - 1; else x2_d = x2 - 1; 
        }
    }
    #ifdef _DEBUG_NEIGHBOURS
    printf("%c-type Neighbours of (%hu, %hu) as determined by dist: \n", synd_type, x1, y1);
    printf("(%hu, %hu), (%hu, %hu)\n", x1_u, y1_u, x1_d, y1_d);
    printf("%c-type Neighbours of (%hu, %hu) as determined by dist: \n", synd_type, x2, y2);
    printf("(%hu, %hu), (%hu, %hu)\n", x2_u, y2_u, x2_d, y2_d);
    #endif
    /*
    Step three: Since the distance between octagonal checks is 
    unambiguous, we determine if any/all of the input co-ordinates 
    represent square checks, and if so, we branch: 
    */
    if (is_sq_cent(x1, y1)){
        #ifdef _DEBUG
        puts("In dist, point 1 is a square");
        #endif
        if (is_sq_cent(x2, y2)){
            #ifdef _DEBUG
            puts("In dist, point 2 is a square");
            printf("Looking through octagonal neighbours: (%hu, %hu), (%hu, %hu), (%hu, %hu), (%hu, %hu) \n",
                (unsigned short int)x1_u, (unsigned short int)y1_u,
                (unsigned short int)x2_u, (unsigned short int)y2_u,
                (unsigned short int)x1_d, (unsigned short int)y1_d,
                (unsigned short int)x2_d, (unsigned short int)y2_d);
            #endif
            num_steps += i_min_4(
                    octagonal_dist(x1_u, y1_u, x2_u, y2_u, sz_x, sz_y),
                    octagonal_dist(x1_u, y1_u, x2_d, y2_d, sz_x, sz_y),
                    octagonal_dist(x1_d, y1_d, x2_u, y2_u, sz_x, sz_y),
                    octagonal_dist(x1_d, y1_d, x2_d, y2_d, sz_x, sz_y)
                    ) + 2;
        }
        else{
            #ifdef _DEBUG
            puts("In dist, point 2 is an octagon");
            #endif
            num_steps += i_min_2(
                        octagonal_dist(x1_u, y1_u, x2, y2, sz_x, sz_y),
                        octagonal_dist(x1_d, y1_d, x2, y2, sz_x, sz_y)
                        ) + 1;
        }
    }
    else if (is_sq_cent(x2, y2)){
        #ifdef _DEBUG
        puts("In dist, point 2 is a square");
        #endif
        num_steps += i_min_2(
                        octagonal_dist(x1, y1, x2_u, y2_u, sz_x, sz_y),
                        octagonal_dist(x1, y1, x2_d, y2_d, sz_x, sz_y)
                        ) + 1;
    }
    else{
        #ifdef _DEBUG
        puts("In dist, both points are octagons");
        #endif
        num_steps += octagonal_dist(x1, y1, x2, y2, sz_x, sz_y);
    }
    return num_steps;
    }
}

num_t toric_dist(num_t x1, num_t y1, num_t x2, num_t y2, num_t sz_x, num_t sz_y, char synd_type){
    /*This function takes the same arguments as squoct_dist, and
    returns the distance between two points on a sz_x by sz_y torus.*/
    return min_abs_sub(x1, x2, sz_x) + min_abs_sub(y1, y2, sz_y);
}