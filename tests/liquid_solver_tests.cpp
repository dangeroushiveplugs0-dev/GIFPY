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
        require(gifpy::LiquidSolver::poly6(0.0, h) > 0.0, "3D poly6 self value must be positive");
        require(near(gifpy::LiquidSolver::poly6(h, h), 0.0), "poly6 must vanish at support radius");
        require(near(gifpy::LiquidSolver::poly6(h + 0.01, h), 0.0), "poly6 must vanish outside support");
        require(near(gifpy::LiquidSolver::viscosity_laplacian(h, h), 0.0), "viscosity kernel support");
        const auto grad = gifpy::LiquidSolver::spiky_gradient({0.0, 0.0, 0.0}, h);
        require(near(grad.x, 0.0) && near(grad.y, 0.0) && near(grad.z, 0.0),
                "gradient at coincident particles must be finite zero");

        gifpy::Parameters p;
        p.gravity = {0.0, -10.0, 0.0};
        p.viscosity = 0.0;
        p.time_step = 0.001;
        p.domain_min = {-100.0, -100.0, -100.0};
        p.domain_max = {100.0, 100.0, 100.0};
        gifpy::LiquidSolver solver(p);
        gifpy::Particle particle;
        particle.position = {0.0, 0.0, 0.0};
        particle.mass = 1.0;
        solver.set_particles({particle});
        solver.step();
        require(near(solver.particles()[0].velocity.y, -0.01, 1.0e-8),
                "single particle gravity update must match semi-implicit Euler");
        require(solver.particles()[0].density > 0.0, "density must include self contribution");

        // Two nearby particles must produce finite, non-identical 3D force updates.
        gifpy::Parameters pair_parameters = p;
        pair_parameters.gravity = {0.0, 0.0, 0.0};
        pair_parameters.rest_density = 1.0;
        pair_parameters.sound_speed = 2.0;
        pair_parameters.smoothing_length = 0.5;
        pair_parameters.time_step = 0.00001;
        gifpy::LiquidSolver pair_solver(pair_parameters);
        gifpy::Particle a;
        a.position = {-0.05, 0.0, 0.0};
        a.velocity = {0.0, 0.0, 0.0};
        a.mass = 0.1;
        gifpy::Particle b = a;
        b.position = {0.05, 0.0, 0.0};
        pair_solver.set_particles({a, b});
        pair_solver.step();
        for (const auto& q : pair_solver.particles()) {
            require(std::isfinite(q.position.x) && std::isfinite(q.position.y) &&
                    std::isfinite(q.position.z) && std::isfinite(q.density) &&
                    std::isfinite(q.pressure), "3D pair interaction must remain finite");
        }
        require(pair_solver.particles()[0].velocity.x < 0.0 &&
                pair_solver.particles()[1].velocity.x > 0.0,
                "positive pressure should separate the particle pair along X");

        // Z-axis wall collision must be handled independently.
        gifpy::Parameters wall_parameters = p;
        wall_parameters.gravity = {0.0, 0.0, 0.0};
        wall_parameters.domain_min = {-1.0, -1.0, 0.0};
        wall_parameters.domain_max = {1.0, 1.0, 1.0};
        wall_parameters.restitution = 1.0;
        gifpy::LiquidSolver wall_solver(wall_parameters);
        gifpy::Particle wall_particle;
        wall_particle.position = {0.0, 0.0, 0.001};
        wall_particle.velocity = {0.0, 0.0, -2.0};
        wall_solver.set_particles({wall_particle});
        wall_solver.step();
        require(wall_solver.particles()[0].position.z >= 0.0,
                "Z boundary must clamp position");
        require(wall_solver.particles()[0].velocity.z > 0.0,
                "Z boundary must reflect outward velocity");

        bool rejected = false;
        try {
            gifpy::Parameters invalid;
            invalid.time_step = 0.0;
            gifpy::LiquidSolver bad(invalid);
        } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "zero timestep must be rejected");

        std::cout << "All GIFPY 3D liquid solver tests passed.\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << "Test failure: " << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
