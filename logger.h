/**
 * @file logger.h
 * @brief The logging subsystem. Only use these functions, the backend may change.
 */
bool log_init(void);
void log_enqueue(int, bool, unsigned long, int);
bool log_framelogger_shutdown(void);
bool log_close(void);
bool log_simulation_stats(unsigned long, double);
void mkpercent(char *, int, double);
