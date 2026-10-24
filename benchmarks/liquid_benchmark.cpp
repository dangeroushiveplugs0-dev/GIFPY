#include "gifpy/liquid_solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::vector<int> parse_counts(const std::string& text) {
    std::vector<int> counts;
    std::stringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        if (item.empty()) throw std::invalid_argument("Empty particle count in --counts");
        std::size_t consumed = 0;
        const long value = std::stol(item, &consumed);
        if (consumed != item.size() || value <= 0 || value > 2000000)
            throw std::invalid_argument("Particle counts must be integers in [1, 2000000]");
        counts.push_back(static_cast<int>(value));
    }
    if (counts.empty()) throw std::invalid_argument("At least one particle count is required");
    std::sort(counts.begin(), counts.end());
    counts.erase(std::unique(counts.begin(), counts.end()), counts.end());
    return counts;
}

std::vector<gifpy::Particle> make_particles(int count, double h) {
    std::vector<gifpy::Particle> particles;
    particles.reserve(static_cast<std::size_t>(count));
    const int side = static_cast<int>(std::ceil(std::cbrt(static_cast<double>(count))));
    for (int i = 0; i < count; ++i) {
        const int x = i % side;
        const int y = (i / side) % side;
        const int z = i / (side * side);
        gifpy::Particle p;
        p.position = {x * h + 0.08 * h * std::sin(i * 1.7),
                      y * h + 0.08 * h * std::sin(i * 0.7 + 1.0),
                      z * h + 0.08 * h * std::cos(i * 1.3)};
        p.mass = 0.01;
        particles.push_back(p);
    }
    return particles;
}
}

int main(int argc, char** argv) {
    try {
        std::vector<int> counts{64, 256, 1024, 4096};
        int steps = 10;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--counts" && i + 1 < argc) {
                counts = parse_counts(argv[++i]);
            } else if (arg == "--steps" && i + 1 < argc) {
                const std::string value = argv[++i];
                std::size_t consumed = 0;
                const long parsed = std::stol(value, &consumed);
                if (consumed != value.size() || parsed <= 0 || parsed > 100000)
                    throw std::invalid_argument("--steps must be an integer in [1, 100000]");
                steps = static_cast<int>(parsed);
            } else {
                throw std::invalid_argument("Usage: gifpy_benchmark [--counts 64,256,1024] [--steps 10]");
            }
        }

        std::cout << "GIFPY 3D WCSPH scaling benchmark\n"
                  << "Deterministic initial state; timings include grid rebuilds and solver steps.\n"
                  << "particles,steps,total_ms,ms_per_step,particle_steps_per_second\n";
        for (const int count : counts) {
            gifpy::Parameters parameters;
            parameters.smoothing_length = 0.1;
            parameters.time_step = 0.0001;
            const double extent = std::max(1.0, std::ceil(std::cbrt(static_cast<double>(count))) * 0.1 + 0.5);
            parameters.domain_min = {-1.0, -1.0, -1.0};
            parameters.domain_max = {extent, extent, extent};
            gifpy::LiquidSolver solver(parameters);
            solver.set_particles(make_particles(count, parameters.smoothing_length));

            const auto start = std::chrono::steady_clock::now();
            for (int step = 0; step < steps; ++step) solver.step();
            const auto finish = std::chrono::steady_clock::now();
            const double elapsed_ms = std::chrono::duration<double, std::milli>(finish - start).count();
            const double step_ms = elapsed_ms / steps;
            const double throughput = elapsed_ms > 0.0
                ? (static_cast<double>(count) * steps * 1000.0) / elapsed_ms : 0.0;
            std::cout << count << ',' << steps << ','
                      << std::fixed << std::setprecision(3) << elapsed_ms << ','
                      << step_ms << ',' << std::setprecision(1) << throughput << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "Benchmark error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
