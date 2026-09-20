#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint>

#include "las-plattice.hpp"
#include "las-reduce-plattice-simd.hpp"

void plattice_info::reduce_many(plattice_info * pli, size_t n, uint32_t I)
{
    size_t i = 0;
#if defined(HAVE_AVX512F)
    for (; i + 16 <= n; i += 16)
        reduce_plattice_simd<16>(pli + i, I);
#endif
#if defined(HAVE_AVX2)
    for (; i + 8 <= n; i += 8)
        reduce_plattice_simd<8>(pli + i, I);
#endif
    for (; i < n; i++)
        pli[i].reduce(I);
}
