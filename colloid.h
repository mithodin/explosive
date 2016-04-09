/**
 ** @file colloid.h
 ** @brief store information for a single colloid
 **/

typedef __m128d vector2d;

/**
 ** @brief Represents one colloid
 ** including position and angle
 **/
typedef struct _colloid {
	vector2d position; /**< x,y position of the colloid */
	double phi; /**< rotation of the colloid */
	struct _colloid *bonding_partner[3];
	int bond_site[3];
	struct _colloid *above;
	struct _colloid *below;
	int external_energy;
} Colloid;

bool colloid_bonded(Colloid *, Colloid *);
