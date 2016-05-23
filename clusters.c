/** @file */
#include <stddef.h>
#include <stdbool.h>
#include <x86intrin.h>
#include "dSFMT/dSFMT.h"
#include "config.h"
#include "geometry.h"
#include "colloid.h"
#include "globals.h"

int iterate_bonding_partners(Colloid *, int, bool*);

/**
 * Calculate the size of the largest cluster
 *
 * @return The size of the largest cluster
 */
int largest_cluster_size(void){
	bool particle_done[NUMBER_OF_PARTICLES]={false};
	int largest_cluster=0;
	int this_cluster=0;
	for(int i=0;i<NUMBER_OF_PARTICLES;++i){
		if( particle_done[i] ) continue;
		this_cluster=iterate_bonding_partners(&(particles[i]),-1,particle_done);
		if(this_cluster > largest_cluster){
			largest_cluster=this_cluster;
		}
		cluster_sizes[this_cluster-1]++;
	}
	return largest_cluster;
}

/**
 * Calculate the size of one cluster recursively
 *
 * @param start The colloid to start iterating the cluster from
 * @param parent The bond site of the colloid I'm coming from (-1 if none)
 * @param particle_done A NUMBER_OF_PARTICLES sized bool array that stores whether or not a particle has already been iterated over
 * @return The size of the cluster
 */
int iterate_bonding_partners(Colloid *start, int parent, bool *particle_done){
	if( particle_done[start->particles_index] ){ return 0; }
	particle_done[start->particles_index]=true;
	int partners=1;
	for(int p=0;p<3;++p){
		if( p != parent && start->bonding_partner[p]!=NULL ){
			partners+=iterate_bonding_partners(start->bonding_partner[p], start->bond_site[p], particle_done);
		}
	}
	return partners;
}
