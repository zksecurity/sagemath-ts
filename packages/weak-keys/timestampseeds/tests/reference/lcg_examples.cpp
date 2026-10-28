// Standalone reference for explicitly named LCG examples, not portable C rand().
// Microsoft CRT constants: UCRT stdlib/rand.cpp; ANSI example: repository README.
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    std::uint32_t state = static_cast<std::uint32_t>(std::strtoull(argv[2], nullptr, 10));
    const auto count = std::strtoul(argv[3], nullptr, 10);
    const bool microsoft = std::string(argv[1]) == "msvc-rand";
    for (unsigned long i = 0; i < count; ++i) {
        state = microsoft ? state * 214013u + 2531011u : (state * 1103515245u + 12345u) & 0x7fffffffu;
        std::cout << (microsoft ? (state >> 16) & 0x7fffu : state) << '\n';
    }
}
