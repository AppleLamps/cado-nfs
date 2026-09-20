#include "cado.h" // IWYU pragma: keep

#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <istream>
#include <string>
#include <vector>

#include "filter_io.hpp"
#include "fstream_maybe_compressed.hpp"
#include "purgedfile.h"
#include "ringbuf.hpp"
#include "timing.h"
#include "typedefs.h"
#include "utils_cxx.hpp"

namespace cado::filter_io_details {
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    static bool filter_rels_force_posix_threads = false;
}
bool cado::filter_io_details::filter_rels_force_posix_threads_is_set()
{
    return filter_rels_force_posix_threads;
}

namespace {
uint32_t read_u32_le(std::istream & f)
{
    unsigned char b[4];
    f.read(reinterpret_cast<char *>(b), 4);
    if (f.gcount() != 4)
        throw cado::error("truncated CADOPURG stream");
    return (uint32_t) b[0]
         | ((uint32_t) b[1] << 8)
         | ((uint32_t) b[2] << 16)
         | ((uint32_t) b[3] << 24);
}

uint64_t read_u64_le(std::istream & f)
{
    uint64_t const lo = read_u32_le(f);
    uint64_t const hi = read_u32_le(f);
    return lo | (hi << 32);
}

int64_t read_i64_le(std::istream & f)
{
    return (int64_t) read_u64_le(f);
}

void append_hex(std::string & s, uint64_t v)
{
    char buf[17];
    snprintf(buf, sizeof(buf), "%" PRIx64, v);
    s += buf;
}

/* Expand a CADOPURG stream (magic already consumed) into the ASCII
 * relation format that the rest of filter_rels parses. */
void feed_purged_binary(ringbuf & r, std::istream & f)
{
    uint32_t const version = read_u32_le(f);
    uint32_t const flags = read_u32_le(f);
    uint64_t const nrows = read_u64_le(f);
    uint64_t const ncols = read_u64_le(f);
    uint64_t const nideals = read_u64_le(f);
    if (version != PURGEDFILE_VERSION)
        throw cado::error("unsupported CADOPURG version {}", version);
#if SIZEOF_INDEX == 8
    if (!(flags & PURGEDFILE_FLAG_INDEX64))
        throw cado::error("CADOPURG file has 32-bit indices; rebuild with SIZEOF_INDEX=4");
#else
    if (flags & PURGEDFILE_FLAG_INDEX64)
        throw cado::error("CADOPURG file has 64-bit indices; rebuild with SIZEOF_INDEX=8");
#endif
    char hdr[128];
    int const nh = snprintf(hdr, sizeof(hdr),
            "# %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",
            nrows, ncols, nideals);
    if (nh > 0)
        r.put(hdr, (size_t) nh);

    std::string line;
    line.reserve(256);
    for (uint64_t i = 0; i < nrows; i++) {
        int64_t const a = read_i64_le(f);
        uint64_t const b = read_u64_le(f);
        uint32_t const n = read_u32_le(f);
        if (n > REL_MAX_SIZE)
            throw cado::error("CADOPURG record with {} ideals (max {})", n, REL_MAX_SIZE);
        line.clear();
        if (a < 0) {
            line.push_back('-');
            append_hex(line, (uint64_t)(-(a + 1)) + 1);
        } else {
            append_hex(line, (uint64_t) a);
        }
        line.push_back(',');
        append_hex(line, b);
        line.push_back(':');
        bool first = true;
        for (uint32_t k = 0; k < n; k++) {
            uint64_t const h = (flags & PURGEDFILE_FLAG_INDEX64)
                ? read_u64_le(f) : read_u32_le(f);
            int c = f.get();
            if (c == EOF)
                throw cado::error("truncated CADOPURG exponent");
            int e = (int8_t) c;
            if (e == 0)
                continue;
            int const times = e < 0 ? -e : e;
            for (int t = 0; t < times; t++) {
                if (!first)
                    line.push_back(',');
                first = false;
                if (e < 0)
                    line.push_back('-');
                append_hex(line, h);
            }
        }
        line.push_back('\n');
        r.put(line.data(), line.size());
    }
}

void feed_one_stream(ringbuf & r, std::istream & f, char const * where)
{
    int const c = f.peek();
    if (c == 'C') {
        char magic[8];
        f.read(magic, sizeof(magic));
        auto const got = f.gcount();
        if (got == 8 && memcmp(magic, PURGEDFILE_MAGIC, 8) == 0) {
            feed_purged_binary(r, f);
            if (!f.good() && !f.eof())
                throw cado::error("read error on {}", where);
            return;
        }
        if (got > 0)
            r.put(magic, (size_t) got);
    }
    r.feed_stream(f);
    if (!f.good() && !f.eof())
        throw cado::error("read error on {}", where);
}
} /* namespace */

void cado::filter_io_details::filter_rels_producer_thread(
    ringbuf & r,
    std::istream& f,
    timingstats_dict_ptr stats)
{
    feed_one_stream(r, f, "input stream");
    r.mark_done();
    if (stats) timingstats_dict_add_mythread(stats, "producer");
}

void cado::filter_io_details::filter_rels_producer_thread(
    ringbuf & r,
    std::vector<std::string> const & input_files,
    timingstats_dict_ptr stats)
{
    for(auto const & filename : input_files) {
        ifstream_maybe_compressed f(filename);
        if (!f)
            throw cado::error("cannot open {}", filename);
        feed_one_stream(r, f, filename.c_str());
    }
    r.mark_done();
    if (stats) timingstats_dict_add_mythread(stats, "producer");
}

void cado::filter_io_details::configure(cxx_param_list & pl)
{
    pl.declare_usage("force-posix-threads", "force the use of posix threads");
    pl.configure_switch("force-posix-threads");
}

void cado::filter_io_details::interpret_parameters(cxx_param_list & pl)
{
    pl.parse("force-posix-threads", filter_rels_force_posix_threads);
}
