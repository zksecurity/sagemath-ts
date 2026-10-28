// CPython integer seed expansion (Modules/_randommodule.c:init_by_array), then
// the platform's std::mt19937 twist/temper engine. No Python runtime is required.
// Seed expansion algorithm: Matsumoto/Nishimura MT19937, as used by CPython.
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>
int main(int argc, char** argv) {
    if (argc != 3) return 2;
    std::vector<std::uint32_t> key(1, 0);
    for (const char digit : std::string(argv[1])) {
        if (digit < '0' || digit > '9') return 2;
        std::uint64_t carry = digit - '0';
        for (auto& word : key) {
            const std::uint64_t value = std::uint64_t(word) * 10 + carry;
            word = std::uint32_t(value); carry = value >> 32;
        }
        if (carry) key.push_back(std::uint32_t(carry));
    }
    std::uint32_t state[624];
    state[0] = 19650218;
    for (std::size_t i = 1; i < 624; ++i) state[i] = 1812433253u * (state[i-1] ^ (state[i-1] >> 30)) + std::uint32_t(i);
    std::size_t i = 1, j = 0;
    for (std::size_t k = key.size() > 624 ? key.size() : 624; k; --k) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1664525u)) + key[j] + std::uint32_t(j);
        ++i; ++j;
        if (i >= 624) { state[0] = state[623]; i = 1; }
        if (j >= key.size()) j = 0;
    }
    for (std::size_t k = 623; k; --k) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1566083941u)) - std::uint32_t(i);
        if (++i >= 624) { state[0] = state[623]; i = 1; }
    }
    state[0] = 0x80000000u;
    std::stringstream serialized;
    for (auto value : state) serialized << value << ' ';
    serialized << 624;
    std::mt19937 engine;
    serialized >> engine;
    if (!serialized) return 2;
    const auto count = std::strtoul(argv[2], nullptr, 10);
    for (unsigned long k = 0; k < count; ++k) std::cout << engine() << '\n';
}
