/**
 * @file logger.h
 * @brief The logging subsystem. Only use these functions, the backend may change.
 */
bool log_init(void);
void log_enqueue(int, bool, unsigned long, int, double);
bool log_framelogger_shutdown(void);
bool log_close(void);
bool log_simulation_stats(unsigned long, double);
void mkpercent(char *, int, double);
bool log_substrate(vector2d *);
bool log_checkpoint(void);
bool log_max_displacement(void);
