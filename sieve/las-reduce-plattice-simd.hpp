#ifndef CADO_LAS_REDUCE_PLATTICE_SIMD_PROD_HPP
#define CADO_LAS_REDUCE_PLATTICE_SIMD_PROD_HPP

/* Production SIMD Franke–Kleinjung reduction. The algorithm is the same
 * as the one exercised in tests/sieve/reduce-plattice/; this header is
 * the copy that las actually calls. */

#include <algorithm>
#include <cstdint>

#include "macros.h"

#include "../tests/sieve/reduce-plattice/simd-base.hpp"
#include "../tests/sieve/reduce-plattice/simd-avx2.hpp"
#include "../tests/sieve/reduce-plattice/simd-avx512.hpp"

template<typename T>
struct simd_div_threshold_prod {
    static constexpr const int value = 4;
};

template<>
struct simd_div_threshold_prod<simd_helper<uint32_t, 8>> {
    static constexpr const int value = 4;
};

#ifdef HAVE_AVX512F
template<>
struct simd_div_threshold_prod<simd_helper<uint32_t, 16>> {
    static constexpr const int value = 5;
};
#endif

template<size_t N>
void reduce_plattice_simd(plattice_info * pli, uint32_t I)
{
    typedef simd_helper<uint32_t, N> A;
    typedef typename A::mask mask;
    typedef typename A::type data;

    data zI = A::set1(I);
    uint32_t explode[N] ATTR_ALIGNED(A::store_alignment);
    for(size_t j = 0 ; j < N ; j++) explode[j] = pli[j].mi0;
    data zmi0 = A::load(explode);
    for(size_t j = 0 ; j < N ; j++) explode[j] = pli[j].i1;
    data zi1 = A::load(explode);
    for(size_t j = 0 ; j < N ; j++) explode[j] = pli[j].j0;
    data zj0 = A::load(explode);
    for(size_t j = 0 ; j < N ; j++) explode[j] = pli[j].j1;
    data zj1 = A::load(explode);

    mask proceed = A::kor(
            A::cmpneq(zj0, A::setzero()),
            A::cmpneq(zj1, A::setzero()));

    mask flip;
    for(flip = A::zeromask() ; ; ) {
        proceed = A::mask_cmpge(proceed, zi1, zI);

        if (!A::mask2int(proceed)) break;

        mask toobig = A::mask_cmpge(proceed, zmi0, A::slli(zi1, simd_div_threshold_prod<A>::value));
        if (UNLIKELY(A::mask2int(toobig))) {
            data k = A::mask_div(A::setzero(), proceed, zmi0, zi1);
            zmi0 = A::sub(zmi0, A::mullo(k, zi1));
            zj0  = A::add(zj0,  A::mullo(k, zj1));
        } else {
            mask subtract = A::mask_cmpge(proceed, zmi0, zi1);
            zmi0 = A::mask_sub(zmi0, subtract, zmi0, zi1);
            zj0  = A::mask_add(zj0,  subtract, zj0,  zj1);
        }
        mask swap = A::cmplt(zmi0, zi1);
        data swapper;
        swapper = A::mask_bxor(A::setzero(), swap, zmi0, zi1);
        zmi0 = A::bxor(zmi0, swapper);
        zi1 = A::bxor(zi1, swapper);
        swapper = A::mask_bxor(A::setzero(), swap, zj0, zj1);
        zj0 = A::bxor(zj0, swapper);
        zj1 = A::bxor(zj1, swapper);
        flip = A::kxor(flip, swap);
    }

    proceed = A::kor(   A::cmpneq(zj0, A::setzero()),
            A::cmpneq(zj1, A::setzero()));

    mask haszero = A::mask_cmpeq(proceed, zi1, A::setzero());

    if (A::mask2int(haszero)) {
        A::store(explode,  zmi0);
        for(size_t j = 0 ; j < N ; j++) pli[j].mi0 = explode[j];
        A::store(explode,  zi1);
        for(size_t j = 0 ; j < N ; j++) pli[j].i1  = explode[j];
        A::store(explode,  zj0);
        for(size_t j = 0 ; j < N ; j++) pli[j].j0  = explode[j];
        A::store(explode,  zj1);
        for(size_t j = 0 ; j < N ; j++) pli[j].j1  = explode[j];

        int p = A::mask2int(proceed);

        for(size_t j = 0, m = A::mask2int(flip) ; j < N ; j++, m>>=1, p>>=1) {
            if (!(p & 1)) continue;
            if (pli[j].i1 == 0) {
                if (!(m&1))
                    pli[j].j0 = pli[j].j1 - pli[j].j0;
                pli[j].reduce_with_vertical_vector(I);
                continue;
            } else {
                int a = (pli[j].mi0 + pli[j].i1 - I) / pli[j].i1;
                pli[j].mi0 -= a * pli[j].i1;
                pli[j].j0  += a * pli[j].j1;
            }
            if (m&1) {
                std::swap(pli[j].mi0, pli[j].i1);
                std::swap(pli[j].j0, pli[j].j1);
            }
        }
    } else {
        data sum = A::sub(A::add(zmi0, zi1), zI);
        mask toobig = A::cmpge(sum, A::slli(zi1, simd_div_threshold_prod<A>::value));
        if (UNLIKELY(A::mask2int(toobig))) {
            data k = A::mask_div(A::setzero(), proceed, sum, zi1);
            zmi0 = A::sub(zmi0, A::mullo(k, zi1));
            zj0  = A::add(zj0,  A::mullo(k, zj1));
        } else {
            for(mask q ; q = A::cmpge(sum, zi1), A::mask2int(q) ; ) {
                zmi0 = A::mask_sub(zmi0, q, zmi0, zi1);
                sum = A::mask_sub(sum, q, sum, zi1);
                zj0 = A::mask_add(zj0, q, zj0, zj1);
            }
        }
        data swapper;
        swapper = A::mask_bxor(A::setzero(), flip, zmi0, zi1);
        zmi0 = A::bxor(zmi0, swapper);
        zi1 = A::bxor(zi1, swapper);
        swapper = A::mask_bxor(A::setzero(), flip, zj0, zj1);
        zj0 = A::bxor(zj0, swapper);
        zj1 = A::bxor(zj1, swapper);

        A::store(explode,  zmi0);
        for(size_t j = 0 ; j < N ; j++) pli[j].mi0 = explode[j];
        A::store(explode,  zi1);
        for(size_t j = 0 ; j < N ; j++) pli[j].i1  = explode[j];
        A::store(explode,  zj0);
        for(size_t j = 0 ; j < N ; j++) pli[j].j0  = explode[j];
        A::store(explode,  zj1);
        for(size_t j = 0 ; j < N ; j++) pli[j].j1  = explode[j];
    }
}

#endif /* CADO_LAS_REDUCE_PLATTICE_SIMD_PROD_HPP */
