/**
 * @file hdf5_logging.h
 * @brief Log information about the simulation in the HDF5 binary format. Do not use directly.
 */
bool h5log_init(void);
bool h5log_log_frame(Colloid *, int, unsigned long, int);
bool h5log_close(void);
bool h5log_log_final_time(unsigned long);
bool h5log_log_acceptance_probability(double);
bool h5log_log_cluster_size(void);
