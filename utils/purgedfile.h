#ifndef CADO_PURGEDFILE_H
#define CADO_PURGEDFILE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Binary purged-file format (little-endian), optional via purge -binary
 * or a .bin / .bin.gz output name. ASCII remains the default.
 *
 *   magic[8]     "CADOPURG"
 *   version      uint32     currently 1
 *   flags        uint32     bit 0: index_t is 64-bit
 *   nrows        uint64     remaining rows
 *   ncols        uint64     column bound (same as ASCII header field 2)
 *   nideals      uint64     remaining columns (ASCII header field 3)
 * then nrows records:
 *   a            int64
 *   b            uint64
 *   n            uint32     number of ideals
 *   h[n]         index_t    (uint32 or uint64 according to flags)
 *   e[n]         int8       exponents
 */
#define PURGEDFILE_MAGIC "CADOPURG"
#define PURGEDFILE_VERSION 1u
#define PURGEDFILE_FLAG_INDEX64 1u

void purgedfile_read_firstline (const char *, uint64_t *, uint64_t *);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_PURGEDFILE_H */
