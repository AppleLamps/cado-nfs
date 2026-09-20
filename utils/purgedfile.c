#include "cado.h" // IWYU pragma: keep

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>      // printf
#include <stdlib.h>        // for abort exit
#include <string.h>

#include "purgedfile.h"
#include "fix-endianness.h"
#include "gzip.h"

/* Read all lines which begin with # in the input file, until we find one
 * which matches the desired format. Binary CADOPURG files are recognized
 * by sniffing the first byte.
 */

// NOLINTBEGIN(clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling,cert-err34-c)
void
purgedfile_read_firstline (const char *fname, uint64_t *nrows, uint64_t *ncols)
{
  FILE *f_tmp = fopen_maybe_compressed (fname, "rb");
  if (!f_tmp)
  {
    fprintf(stderr, "%s: %s\n", fname, strerror(errno));
    abort();
  }
  int c = fgetc (f_tmp);
  if (c == EOF) {
    fprintf(stderr,
            "Parse error while reading %s: empty file\n", fname);
    fclose_maybe_compressed(f_tmp, fname);
    exit(EXIT_FAILURE);
  }
  ungetc (c, f_tmp);
  if (c == 'C') {
    char magic[8];
    if (fread (magic, 1, 8, f_tmp) != 8
        || memcmp (magic, PURGEDFILE_MAGIC, 8) != 0) {
      fprintf(stderr,
              "Parse error while reading %s: not a CADOPURG file\n", fname);
      fclose_maybe_compressed(f_tmp, fname);
      exit(EXIT_FAILURE);
    }
    uint32_t version = 0, flags = 0;
    uint64_t nideals = 0;
    if (fread32_little (&version, 1, f_tmp) != 1
        || fread32_little (&flags, 1, f_tmp) != 1
        || fread64_little (nrows, 1, f_tmp) != 1
        || fread64_little (ncols, 1, f_tmp) != 1
        || fread64_little (&nideals, 1, f_tmp) != 1) {
      fprintf(stderr,
              "Parse error while reading %s: truncated CADOPURG header\n", fname);
      fclose_maybe_compressed(f_tmp, fname);
      exit(EXIT_FAILURE);
    }
    if (version != PURGEDFILE_VERSION) {
      fprintf(stderr,
              "Parse error while reading %s: unsupported CADOPURG version %" PRIu32 "\n",
              fname, version);
      fclose_maybe_compressed(f_tmp, fname);
      exit(EXIT_FAILURE);
    }
    (void) flags;
    (void) nideals;
    fclose_maybe_compressed(f_tmp, fname);
    return;
  }
  char buf[1024];
  while (fgets(buf, sizeof(buf), f_tmp)) {
      if (*buf != '#') {
          break;
      }
      int ret = sscanf(buf, "# %" SCNu64 " %" SCNu64 "", nrows, ncols);
      if (ret == 2) {
          fclose_maybe_compressed(f_tmp, fname);
          return;
      }
  }
  fprintf(stderr,
          "Parse error while reading %s: no header line with desired format\n", fname);
  fclose_maybe_compressed(f_tmp, fname);
  exit(EXIT_FAILURE);
}
// NOLINTEND(clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling,cert-err34-c)
