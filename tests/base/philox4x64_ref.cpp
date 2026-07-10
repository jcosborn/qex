#include <Random123/philox.h>
#include <Random123/uniform.hpp>

#include <cerrno>
#include <cctype>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

// clang++ -std=c++11 -Iextra/Random123/include \
//   tests/base/philox4x64_ref.cpp -o philox4x64_ref

static_assert(philox4x64_rounds == 10,
              "the reference output requires Philox4x64/10");
static_assert(sizeof(double) == sizeof(uint64_t) &&
              std::numeric_limits<double>::is_iec559 &&
              std::numeric_limits<double>::digits == 53,
              "the reference output requires IEEE binary64");

struct Engine {
    philox4x64_ctr_t ctr;
    philox4x64_key_t key;
    philox4x64_ctr_t out;
    unsigned int part;

    Engine(philox4x64_ctr_t c, philox4x64_key_t k, unsigned int p)
        : ctr(c), key(k), out(philox4x64(ctr, key)), part(p) {}

    uint32_t next() {
        const uint64_t x = out.v[part >> 1];
        const uint32_t r = (part & 1) ? uint32_t(x >> 32) : uint32_t(x);
        if (++part == 8) {
            part = 0;
            ctr.incr();
            out = philox4x64(ctr, key);
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
    if (s[0] == '-' || std::isspace(static_cast<unsigned char>(s[0])))
        return false;
    char *end;
    errno = 0;
    const unsigned long long y = std::strtoull(s, &end, 0);
    if (errno != 0 || end == s || *end != '\0' ||
        y > std::numeric_limits<uint64_t>::max())
        return false;
    x = uint64_t(y);
    return true;
}

int main(int argc, char **argv) {
    if (argc < 7 || argc > 9) {
        std::fprintf(stderr,
            "usage: %s c0 c1 c2 c3 k0 k1 [part [count]]\n", argv[0]);
        return 1;
    }

    uint64_t a[8] = {};
    for (int i = 1; i < argc; ++i) {
        if (!arg(argv[i], a[i-1])) {
            std::fprintf(stderr, "invalid argument: %s\n", argv[i]);
            return 1;
        }
    }
    philox4x64_ctr_t ctr = {{a[0], a[1], a[2], a[3]}};
    philox4x64_key_t key = {{a[4], a[5]}};
    const uint64_t part64 = argc > 7 ? a[6] : 0;
    const uint64_t count = argc > 8 ? a[7] : 4;
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
