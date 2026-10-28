// Independent oracle: the platform C++ standard library, integer seeding only.
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    std::mt19937 generator(static_cast<std::uint32_t>(std::strtoull(argv[1], nullptr, 10)));
    const auto count = std::strtoul(argv[2], nullptr, 10);
    for (unsigned long i = 0; i < count; ++i) std::cout << generator() << '\n';
}
