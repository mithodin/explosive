/** @file */

#include <string.h>
#include <limits.h>     /* PATH_MAX */
#include <sys/stat.h>   /* mkdir(2) */
#include <errno.h>
#include "helpers.h"

/**
 * Create path. Equivalent to mkdir -p $path. New directories will be created mode 770
 * @param path The path to create
 * @return 0 on success, -1 otherwise
 */
int mkdir_p(const char *path)
{
	const size_t len = strlen(path);
	char _path[PATH_MAX];
	char *p; 

	errno = 0;

	/* Copy string so its mutable */
	if (len > sizeof(_path)-1) {
		errno = ENAMETOOLONG;
		return -1; 
	}   
	strcpy(_path, path);

	/* Iterate the string */
	for (p = _path + 1; *p; p++) {
		if (*p == '/') {
			/* Temporarily truncate */
			*p = '\0';

			if (mkdir(_path, S_IRWXU | S_IRWXG) != 0) {
				if (errno != EEXIST)
					return -1; 
			}

			*p = '/';
		}
	}   

	if (mkdir(_path, S_IRWXU) != 0) {
		if (errno != EEXIST)
			return -1; 
	}   

	return 0;
}
