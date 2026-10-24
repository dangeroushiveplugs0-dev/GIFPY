#include "gifpy/liquid_solver.hpp"

#include <iostream>
#include <vector>

int main() {
    gifpy::Parameters parameters;
    parameters.smoothing_length = 0.15;
    parameters.time_step = 0.0001;
    parameters.domain_min = {0.0, 0.0};
    parameters.domain_max = {1.0, 1.0};

    gifpy::LiquidSolver solver(parameters);
    std::vector<gifpy::Particle> particles;
    for (int y = 0; y < 5; ++y) {
        for (int x = 0; x < 5; ++x) {
            gifpy::Particle p;
            p.position = {0.35 + 0.04 * x, 0.55 + 0.04 * y};
            p.mass = 0.01;
            particles.push_back(p);
        }
    }
    solver.set_particles(std::move(particles));

    for (int i = 0; i < 100; ++i) solver.step();
    const auto& result = solver.particles();
    std::cout << "GIFPY WCSPH smoke test: " << result.size()
              << " particles after 100 steps\n";
    for (std::size_t i = 0; i < result.size(); i += 6) {
        std::cout << i << ": x=(" << result[i].position.x << ", "
                  << result[i].position.y << "), rho=" << result[i].density << '\n';
    }
    return 0;
}
