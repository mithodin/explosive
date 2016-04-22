/**
 * @file colloid.h
 * @brief Data types and functions concerning single colloids
 **/

/**
 * @brief Represents one colloid
 * including position and angle
 **/
typedef struct _colloid {
	vector2d position; /**< x,y position of the colloid */
	double phi; /**< rotation of the colloid */
	struct _colloid *bonding_partner[3]; /**< current bonding partners */
	int bond_site[3]; /**< which patch on the partner am I bonded to? */
	struct _colloid *above; /**< which colloid is above me? */
	struct _colloid *below; /**< which colloid is below me? */
	double external_energy; /**< current external energy (substrate) */
	int internal_energy; /**< current internal energy (bonds) */
	int particles_index; /**< where is this colloid in the particles array? */
} Colloid;

bool colloid_bonded(Colloid *, Colloid *, bool*, int*, int*);
void init_ysorted_list(void);
void insert_sorted_y(Colloid *, Colloid *);
void init_bonding_partners(void);

/** initializer for an empty colloid */
#define EMPTY_COLLOID {{0.0,0.0},0.0,{NULL,NULL,NULL},{-1,-1,-1},NULL,NULL,0,0,0}
