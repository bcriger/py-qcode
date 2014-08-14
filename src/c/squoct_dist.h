#ifndef squoct_dist_h__
#define squoct_dist_h__

typedef unsigned short int num_t;
typedef enum { false, true } bool;

extern num_t dist(num_t x1, num_t y1, num_t x2, num_t y2, num_t sz_x, num_t sz_y, char synd_type);

#endif