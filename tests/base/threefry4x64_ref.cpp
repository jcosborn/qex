#include <Random123/threefry.h>
#include <Random123/uniform.hpp>

#include <cerrno>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

// clang++ -std=c++11 -Iextra/Random123/include \
//   tests/base/threefry4x64_ref.cpp -o threefry4x64_ref

static_assert(sizeof(double) == sizeof(uint64_t) &&
              std::numeric_limits<double>::is_iec559 &&
              std::numeric_limits<double>::digits == 53,
              "the reference output requires IEEE binary64");

struct Engine {
    threefry4x64_ctr_t ctr;
    threefry4x64_key_t key;
    threefry4x64_ctr_t out;
    unsigned int part;

    Engine(threefry4x64_ctr_t c, threefry4x64_key_t k, unsigned int p)
        : ctr(c), key(k), out(threefry4x64(ctr, key)), part(p) {}

    uint32_t next() {
        const uint64_t x = out.v[part >> 1];
        const uint32_t r = (part & 1) ? uint32_t(x >> 32) : uint32_t(x);
        if (++part == 8) {
            part = 0;
            ctr.incr();
            out = threefry4x64(ctr, key);
        }
        return r;
    }

    uint64_t next64() {
        const uint64_t lo = next();
        const uint64_t hi = next();
        return lo | (hi << 32);
    }
};

static bool arg(const char *s, uint64_t& x) {
    char *end;
    errno = 0;
    const unsigned long long y = std::strtoull(s, &end, 0);
    if (s[0] == '-' || errno != 0 || end == s || *end != '\0' ||
        y > std::numeric_limits<uint64_t>::max())
        return false;
    x = uint64_t(y);
    return true;
}

int main(int argc, char **argv) {
    if (argc < 9 || argc > 11) {
        std::fprintf(stderr,
            "usage: %s c0 c1 c2 c3 k0 k1 k2 k3 [part [count]]\n", argv[0]);
        return 1;
    }

    uint64_t a[10] = {};
    for (int i = 1; i < argc; ++i) {
        if (!arg(argv[i], a[i-1])) {
            std::fprintf(stderr, "invalid argument: %s\n", argv[i]);
            return 1;
        }
    }
    threefry4x64_ctr_t ctr = {{a[0], a[1], a[2], a[3]}};
    threefry4x64_key_t key = {{a[4], a[5], a[6], a[7]}};
    const uint64_t part64 = argc > 9 ? a[8] : 0;
    const uint64_t count = argc > 10 ? a[9] : 4;
    if (part64 >= 8) {
        std::fprintf(stderr, "part must be in 0..7\n");
        return 1;
    }

    Engine rng(ctr, key, unsigned(part64));
    for (uint64_t i = 0; i < count; ++i) {
        const uint64_t x = rng.next64();
        const double u = r123::u01<double>(x);
        uint64_t bits;
        std::memcpy(&bits, &u, sizeof(bits));
        std::printf("%016" PRIx64 " %a %016" PRIx64 "\n", x, u, bits);
    }
    return 0;
}
