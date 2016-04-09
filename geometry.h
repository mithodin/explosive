/**
 * @file geometry.h
 * @brief Collection of functions dealing with the geometry of the box
 */

double distance_direct(double, double, double, double);
vector2d distance(vector2d, vector2d, double*);
void make_periodic(vector2d *);
double angle_twopi(double);
