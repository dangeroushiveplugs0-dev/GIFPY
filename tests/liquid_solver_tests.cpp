#include "gifpy/liquid_solver.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
bool near(double a, double b, double tolerance = 1.0e-9) {
    return std::abs(a - b) <= tolerance;
}
}

int main() {
    try {
        const double h = 0.2;
        require(gifpy::LiquidSolver::poly6(0.0, h) > 0.0, "poly6 self value must be positive");
        require(near(gifpy::LiquidSolver::poly6(h, h), 0.0), "poly6 must vanish at support radius");
        require(near(gifpy::LiquidSolver::poly6(h + 0.01, h), 0.0), "poly6 must vanish outside support");
        require(near(gifpy::LiquidSolver::viscosity_laplacian(h, h), 0.0), "viscosity kernel support");
        const auto grad = gifpy::LiquidSolver::spiky_gradient({0.0, 0.0}, h);
        require(near(grad.x, 0.0) && near(grad.y, 0.0), "gradient at coincident particles must be finite zero");

        gifpy::Parameters p;
        p.gravity = {0.0, -10.0};
        p.viscosity = 0.0;
        p.time_step = 0.001;
        p.domain_min = {-100.0, -100.0};
        p.domain_max = {100.0, 100.0};
        gifpy::LiquidSolver solver(p);
        gifpy::Particle particle;
        particle.position = {0.0, 0.0};
        particle.mass = 1.0;
        solver.set_particles({particle});
        solver.step();
        require(near(solver.particles()[0].velocity.y, -0.01, 1.0e-8),
                "single particle gravity update must match semi-implicit Euler");
        require(solver.particles()[0].density > 0.0, "density must include self contribution");

        bool rejected = false;
        try {
            gifpy::Parameters invalid;
            invalid.time_step = 0.0;
            gifpy::LiquidSolver bad(invalid);
        } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "zero timestep must be rejected");

        std::cout << "All GIFPY liquid solver tests passed.\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << "Test failure: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
