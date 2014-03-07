typedef enum {TWOPATCH,THREEPATCH} species;

struct colloid {
	double x; //Coordinates
	double z;
	double a;
	species sp; //what species?
	struct colloid *left; //Keeps a sorted double linked list in x-direction
	struct colloid *right;
	struct colloid *above; //and one in z-direction
	struct colloid *below;
};

typedef struct colloid Colloid;

void insertSortedX(Colloid *list, Colloid *newitem);
void insertSortedZ(Colloid *list, Colloid *newitem);
void insertLeft(Colloid *list, Colloid *newitem);
void insertRight(Colloid *list, Colloid *newitem);
void newColloid(species, Colloid *);
void printColloidsSortedX(Colloid *);
