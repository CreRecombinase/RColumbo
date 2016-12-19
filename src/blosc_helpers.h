#ifndef BLOSC_HELPERS_H
#define BLOSC_HELPERS_H
#include "blosc_filter.h"
#include <hdf5.h>

#if defined(__GNUC__)
#define PUSH_ERR(func, minor, str, ...) H5Epush(H5E_DEFAULT, __FILE__, func, __LINE__, H5E_ERR_CLS, H5E_PLINE, minor, str, ##__VA_ARGS__)
#elif defined(_MSC_VER)
#define PUSH_ERR(func, minor, str, ...) H5Epush(H5E_DEFAULT, __FILE__, func, __LINE__, H5E_ERR_CLS, H5E_PLINE, minor, str, __VA_ARGS__)
#else
/* This version is portable but it's better to use compiler-supported
 approaches for handling the trailing comma issue when possible. */
#define PUSH_ERR(func, minor, ...) H5Epush(H5E_DEFAULT, __FILE__, func, __LINE__, H5E_ERR_CLS, H5E_PLINE, minor, __VA_ARGS__)
#endif	/* defined(__GNUC__) */

#define GET_FILTER(a,b,c,d,e,f,g) H5Pget_filter_by_id(a,b,c,d,e,f,g,NULL)


size_t blosc_filter(unsigned flags, size_t cd_nelmts,
                    const unsigned cd_values[], size_t nbytes,
                    size_t *buf_size, void **buf);

herr_t blosc_set_local(hid_t dcpl, hid_t type, hid_t space);


/* Register the filter, passing on the HDF5 return value */

/*  Filter setup.  Records the following inside the DCPL:

 1. If version information is not present, set slots 0 and 1 to the filter
 revision and Blosc version, respectively.

 2. Compute the type size in bytes and store it in slot 2.

 3. Compute the chunk size in bytes and store it in slot 3.
 */



#endif
