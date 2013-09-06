/*
 * sys/mman.h
 * mman-win32
 */


//#ifndef _SYS_MMAN_H_
//#define _SYS_MMAN_H_
#ifndef _MMAP_WIN_H
#define _MMAP_WIN_H

#include <stdio.h>
#include <windows.h>
#include <sys/types.h>


#ifdef __cplusplus
extern "C" {
#endif



//#define PROT_READ       1
//#define PROT_WRITE      2
//#define PROT_EXEC       4


//#define MAP_FILE        0
//#define MAP_SHARED      1
//#define MAP_PRIVATE     2
//#define MAP_TYPE        0xf
#define MAP_FIXED         0x10
//#define MAP_ANONYMOUS   0x20
//#define MAP_ANON        MAP_ANONYMOUS


#ifdef _SYS_MMAN_H_
#undef PROT_READ
#undef PROT_WRITE
#undef PROT_EXEC
#undef MAP_SHARED
#undef MAP_PRIVATE
#undef MAP_FAILED
#endif

#define PROT_NONE        0x0
#define PROT_READ        0x1
#define PROT_WRITE       0x2

/* This flag is only available in WinXP+ */
#ifdef FILE_MAP_EXECUTE
#define PROT_EXEC        0x4
#else
#define PROT_EXEC        0x0
#define FILE_MAP_EXECUTE 0
#endif

#define MAP_SHARED       0x01
#define MAP_PRIVATE      0x02

/* Prefer MAP_ANONYMOUS since MAP_ANON is deprecated according to man page. */
#if !defined(MAP_ANONYMOUS) && defined(MAP_ANON)
#  define MAP_ANONYMOUS MAP_ANON
#else
#  define MAP_ANONYMOUS 0x20
#  define MAP_ANON      MAP_ANONYMOUS
#endif


#define MAP_FAILED      ((void *)-1)


/* Flags for msync. */
#define MS_ASYNC        1
#define MS_SYNC         2
#define MS_INVALIDATE   4


void*   mmap(void* addr, size_t len, int prot, int flags, int fildes, off_t off);
int     munmap(void* addr, size_t len);
int     mprotect(void* addr, size_t len, int prot);
int     msync(void* addr, size_t len, int flags);
int     mlock(const void* addr, size_t len);
int     munlock(const void* addr, size_t len);


#ifdef __cplusplus
};
#endif


#endif /*  _SYS_MMAN_H_ */