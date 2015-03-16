"""
This script will call the function 'dist' from libqcode_dist.so in
../c, using ctypes.
"""
from ctypes import CDLL, c_ushort, c_char

libqcode_dist = CDLL('../c/libqcode_dist.so')

total_size = (24, 24)

# Convenience function for checking distances:


def print_dist(pt1, pt2):
    for synd_type in 'XZ':
        x1, y1 = map(c_ushort, pt1)
        x2, y2 = map(c_ushort, pt2)
        sz_x, sz_y = map(c_ushort, total_size)
        delta = libqcode_dist.squoct_dist(x1, y1, x2, y2, sz_x, sz_y,
                                          c_char(synd_type))
        print "{3}-type distance between {0} and {1}: {2}".format(
            pt1, pt2, delta, synd_type)
    pass

"""
Initialize a set of test point pairs, each pair tests a different
aspect of the function.
"""

# 2D Octagon-Octagon wrapping in both directions:
pt1 = (4, 4)
pt2 = (22, 22)
print_dist(pt1, pt2)

# 2D Square-Octagon wrapping in both directions:
pt1 = (7, 4)
pt2 = (22, 22)
print_dist(pt1, pt2)

# 2D Square-Square wrapping in both directions:
pt1 = (7, 4)
pt2 = (7, 13)
print_dist(pt1, pt2)
