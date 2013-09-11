#ifndef _UNISTD_H
#define _UNISTD_H        1

/* This file intended to serve as a drop-in replacement for
 *  *  unistd.h on Windows
 *   *  Please add functionality as neeeded
 *    */

//#include <stdlib.h>
//#include <io.h>

//#define srandom srand
//#define random rand

/* This implementation is only for native Win32 systems.  */
#if (defined _WIN32 || defined __WIN32__) && ! defined __CYGWIN__

# define WIN32_LEAN_AND_MEAN
# include <windows.h>

int
getpagesize(void)
{
    SYSTEM_INFO system_info;
    GetSystemInfo(&system_info);
    return system_info.dwPageSize;
}

#endif


//#define access _access
#define ftruncate _chsize

//#define ssize_t int

#endif /* unistd.h  */
