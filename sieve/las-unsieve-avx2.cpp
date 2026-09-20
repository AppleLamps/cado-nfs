#include "cado.h" // IWYU pragma: keep

#ifndef HAVE_AVX2
#error "This file assumes AVX2 support!"
#endif /* HAVE_AVX2 */

#include <cstdint>
#include <cstdlib>

#include <vector>

#include <immintrin.h>

#ifdef TRACE_K
#include "las-where-am-i.hpp"
#include "las-output.hpp"
#include "verbose.hpp"
#endif

#include "las-unsieve.hpp"
#include "macros.h"
#include "arithxx/u64arith.h"
#include "gcd.h"

static const int verify_gcd = 0;
alignas(32) static const __m256i sign_conversion = _mm256_set1_epi8(-128);
alignas(32) static const __m256i ff = _mm256_set1_epi8(0xff);
/* Byte 0 is the lowest address. When j is even, even i (even x if i0 is
 * even) must not survive: byte 0,2,4,... are 0. Same layout as the SSE2
 * even_masks, just 32 bytes wide. */
alignas(32) static const __m256i even_masks[2] = {
    _mm256_set_epi8(0xFF, 0, 0xFF, 0, 0xFF, 0, 0xFF, 0, 0xFF, 0, 0xFF, 0,
                    0xFF, 0, 0xFF, 0, 0xFF, 0, 0xFF, 0, 0xFF, 0, 0xFF, 0,
                    0xFF, 0, 0xFF, 0, 0xFF, 0, 0xFF, 0),
    _mm256_set1_epi8(0xff)};

static inline uint64_t
sieve_info_test_lognorm_avx2_mask(unsigned char * S0, const __m256i pattern0,
                                  unsigned char const * S1, const __m256i pattern1)
{
    __m256i const a = _mm256_loadu_si256((__m256i const *) S0);
    __m256i const r = _mm256_loadu_si256((__m256i const *) S1);
    __m256i m1 = _mm256_cmpgt_epi8(pattern0, _mm256_xor_si256(a, sign_conversion));
    __m256i const m2 = _mm256_cmpgt_epi8(pattern1, _mm256_xor_si256(r, sign_conversion));
    m1 = _mm256_and_si256(m1, m2);
    _mm256_storeu_si256((__m256i *) S0, _mm256_or_si256(a, _mm256_xor_si256(m1, ff)));
    /* movemask is a 32-bit int. Zero-extend so that a survivor in byte 31
     * can be shifted out without a 32-bit shift-by-32 (undefined). */
    return (uint64_t) (unsigned int) _mm256_movemask_epi8(m1);
}

static inline uint64_t
sieve_info_test_lognorm_avx2_mask_oneside(unsigned char * S0, const __m256i pattern0)
{
    __m256i const a = _mm256_loadu_si256((__m256i const *) S0);
    __m256i const m1 = _mm256_cmpgt_epi8(pattern0, _mm256_xor_si256(a, sign_conversion));
    _mm256_storeu_si256((__m256i *) S0, _mm256_or_si256(a, _mm256_xor_si256(m1, ff)));
    return (uint64_t) (unsigned int) _mm256_movemask_epi8(m1);
}

static inline void
search_single_survivors_mask_avx2(unsigned char * const SS,
        unsigned int j,
        int i0,
        int i1 MAYBE_UNUSED,
        int N MAYBE_UNUSED,
        int x_start,
        unsigned int nr_div,
        unsigned int (*div)[2],
        uint64_t bitmask,
        std::vector<uint32_t> &survivors)
{
    for (int x = x_start; UNLIKELY(bitmask != 0); x++) {
        const unsigned int tz = u64arith_ctz(bitmask);
        x += (int) tz;
        bitmask >>= tz + 1u;

        const unsigned int i = abs(i0 + x);
        int divides = 0;
        switch (nr_div) {
            // coverity[unterminated_case]
          case 6: divides |= (i * div[5][0] <= div[5][1]); no_break();
            // coverity[unterminated_case]
          case 5: divides |= (i * div[4][0] <= div[4][1]); no_break();
            // coverity[unterminated_case]
          case 4: divides |= (i * div[3][0] <= div[3][1]); no_break();
            // coverity[unterminated_case]
          case 3: divides |= (i * div[2][0] <= div[2][1]); no_break();
            // coverity[unterminated_case]
          case 2: divides |= (i * div[1][0] <= div[1][1]); no_break();
            // coverity[unterminated_case]
          case 1: divides |= (i * div[0][0] <= div[0][1]); no_break();
          case 0: while(0){};
        }

        if (divides)
        {
            if (verify_gcd)
                ASSERT_ALWAYS(bin_gcd_int64_safe (i, j) != 1);
#ifdef TRACE_K
            if (trace_on_spot_Nx(N, x)) {
                verbose_fmt_print(TRACE_CHANNEL, 0, "# Slot [{}] in bucket {} has non coprime (i,j)=({},{})\n",
                        x, N, i, j);
            }
#endif
            SS[x] = 255;
        } else {
            survivors.push_back(x);
            if (verify_gcd)
                ASSERT_ALWAYS(bin_gcd_int64_safe (i, j) == 1);
#ifdef TRACE_K
            if (trace_on_spot_Nx(N, x)) {
                verbose_fmt_print(TRACE_CHANNEL, 0, "# Slot [{}] in bucket {} is survivor with coprime (i,j)\n",
                        x, N);
            }
#endif
        }
    }
}

static void
search_survivors_in_line1_avx2(unsigned char * const SS[2],
        const unsigned char bound[2],
        unsigned int j,
        int i0, int i1,
        int N,
        j_divisibility_helper const & j_div,
        unsigned int td_max,
        std::vector<uint32_t> &survivors)
{
    unsigned int div[6][2], nr_div;
    nr_div = extract_j_div(div, j, j_div, 3, td_max);
    ASSERT_ALWAYS(nr_div <= 6);

    const __m256i even_mask = even_masks[j % 2];
    __m256i const patterns[2] = {
        _mm256_xor_si256(_mm256_and_si256(_mm256_set1_epi8(bound[0] + 1), even_mask), sign_conversion),
        _mm256_xor_si256(_mm256_and_si256(_mm256_set1_epi8(bound[1] + 1), even_mask), sign_conversion)
    };
    const int x_step = (int) sizeof(__m256i);

    for (int x_start = 0; x_start < (i1 - i0); x_start += x_step)
    {
        const uint64_t mask = sieve_info_test_lognorm_avx2_mask(
                    SS[0] + x_start, patterns[0],
                    SS[1] + x_start, patterns[1]);
        search_single_survivors_mask_avx2(SS[0], j, i0, i1, N, x_start,
            nr_div, div, mask, survivors);
    }
}

static void
search_survivors_in_line1_avx2_oneside(unsigned char * SS,
        const unsigned char bound,
        unsigned int j,
        int i0, int i1,
        int N,
        j_divisibility_helper const & j_div,
        unsigned int td_max,
        std::vector<uint32_t> &survivors)
{
    unsigned int div[6][2], nr_div;
    nr_div = extract_j_div(div, j, j_div, 3, td_max);
    ASSERT_ALWAYS(nr_div <= 6);

    const __m256i even_mask = even_masks[j % 2];
    __m256i const pattern =
        _mm256_xor_si256(_mm256_and_si256(_mm256_set1_epi8(bound + 1), even_mask), sign_conversion);
    const int x_step = (int) sizeof(__m256i);

    for (int x_start = 0; x_start < (i1 - i0); x_start += x_step)
    {
        const uint64_t mask = sieve_info_test_lognorm_avx2_mask_oneside(
                    SS + x_start, pattern);
        search_single_survivors_mask_avx2(SS, j, i0, i1, N, x_start,
            nr_div, div, mask, survivors);
    }
}

void
search_survivors_in_line_avx2(unsigned char * const SS[2],
        const unsigned char bound[2],
        unsigned int j,
        int i0, int i1,
        int N,
        j_divisibility_helper const & j_div,
        unsigned int td_max, std::vector<uint32_t> &survivors)
{
    /* Pattern-sieve of 3|gcd(i,j) and 5|gcd(i,j) still uses the SSE2
     * 16-byte patterns stored in unsieve_data. The common path (most
     * lines) is AVX2. */
#ifdef HAVE_SSE2
    if (j % 3 == 0 || j % 5 == 0 || ((i1 - i0) & 31) != 0) {
        search_survivors_in_line_sse2(SS, bound, j, i0, i1, N, j_div,
                td_max, survivors);
        return;
    }
#endif
    search_survivors_in_line1_avx2(SS, bound, j, i0, i1, N, j_div,
            td_max, survivors);
}

void
search_survivors_in_line_avx2_oneside(unsigned char * const SS,
        const unsigned char bound,
        unsigned int j,
        int i0, int i1,
        int N,
        j_divisibility_helper const & j_div,
        unsigned int td_max, std::vector<uint32_t> &survivors)
{
#ifdef HAVE_SSE2
    if (j % 3 == 0 || j % 5 == 0 || ((i1 - i0) & 31) != 0) {
        search_survivors_in_line_sse2_oneside(SS, bound, j, i0, i1, N, j_div,
                td_max, survivors);
        return;
    }
#endif
    search_survivors_in_line1_avx2_oneside(SS, bound, j, i0, i1, N, j_div,
            td_max, survivors);
}

void
search_survivors_in_line_avx2_siqs(
        unsigned char * SS,
        unsigned char bound,
        unsigned int length,
        std::vector<uint16_t> &survivors)
{
    __m256i const B = _mm256_xor_si256(_mm256_set1_epi8(bound+1), sign_conversion);
    const unsigned int x_step = sizeof(__m256i);

#ifdef HAVE_SSE2
    if ((length & 31u) != 0) {
        search_survivors_in_line_sse2_siqs(SS, bound, length, survivors);
        return;
    }
#endif
    for (unsigned int x_start = 0; x_start < length; x_start += x_step)
    {
        /* Do bounds check using AVX2 pattern, set non-survivors in SS[0] array
           to 255 */
        __m256i const s = _mm256_loadu_si256((__m256i const *)(SS + x_start));
        __m256i m = _mm256_cmpgt_epi8(B, _mm256_xor_si256(s, sign_conversion));
        /* movemask returns an int, we want to zero-extend the result into an
         * uint64_t. To avoid sign-extension we need to first cast it into an
         * unsigned int and then into an uint64_t */
        uint64_t bitmask = (uint64_t) (unsigned int) _mm256_movemask_epi8(m);
        m = _mm256_xor_si256(m, ff);
        _mm256_storeu_si256((__m256i *)(SS + x_start), _mm256_or_si256(s, m));

        for (unsigned int x = x_start; UNLIKELY(bitmask != 0); ++x) {
            unsigned int const tz = u64arith_ctz(bitmask);
            x += tz;
            bitmask >>= tz + 1u;

            survivors.push_back(x);
        }
    }
}
