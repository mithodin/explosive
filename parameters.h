#define SIZE_X 2.0
#define SIZE_Y 2.0
#define PERIODIC_X
#define PERIODIC_Y

extern const double colloid_diameter; /**< Diameter of a colloid. Defines the unit of length */
extern const double colloid_patch_diameter; /**< Diameter of a patch on the surface of the colloid */
extern const double colloid_min_bond_dist; /**< Minimum distance between colloids where a bond is possible */
extern const double energy_well_depth;
extern const double energy_bond;
extern const int num_particles;

extern const double onethirdpi;
extern const double twothirdspi;

extern dsfmt_t rng;
