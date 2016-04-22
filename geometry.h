/**
 * @file geometry.h
 * @brief Collection of functions dealing with the geometry of the box
 */

/**
 * @brief two double-precision floats that can be worked with using SIMD instructions
 */
typedef __m128d vector2d;
typedef __m256d vector4d;

double distance_direct(double, double, double, double);
vector2d distance(vector2d, vector2d, double*);
void make_periodic(vector2d *);
double angle_twopi(double);
double distance_y(double, double);
