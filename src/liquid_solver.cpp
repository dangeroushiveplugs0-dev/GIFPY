#include "gifpy/liquid_solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace gifpy {
namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double density_floor = 1.0e-9;

bool finite(Vec3 v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}
}

Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
Vec3 operator-(Vec3 a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
Vec3 operator*(Vec3 a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}
double length(Vec3 v) { return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z); }

LiquidSolver::LiquidSolver(Parameters parameters) : parameters_(parameters) {
    validate_parameters();
}

void LiquidSolver::validate_parameters() const {
    const auto finite_scalar = [](double v) { return std::isfinite(v); };
    if (!finite_scalar(parameters_.smoothing_length) || parameters_.smoothing_length <= 0.0 ||
        !finite_scalar(parameters_.rest_density) || parameters_.rest_density <= 0.0 ||
        !finite_scalar(parameters_.sound_speed) || parameters_.sound_speed <= 0.0 ||
        !finite_scalar(parameters_.viscosity) || parameters_.viscosity < 0.0 ||
        !finite_scalar(parameters_.time_step) || parameters_.time_step <= 0.0 ||
        !finite_scalar(parameters_.restitution) || parameters_.restitution < 0.0 ||
        parameters_.restitution > 1.0 ||
        !finite(parameters_.gravity) || !finite(parameters_.domain_min) ||
        !finite(parameters_.domain_max) ||
        parameters_.domain_min.x >= parameters_.domain_max.x ||
        parameters_.domain_min.y >= parameters_.domain_max.y ||
        parameters_.domain_min.z >= parameters_.domain_max.z) {
        throw std::invalid_argument("Invalid GIFPY 3D liquid solver parameters");
    }
}

void LiquidSolver::set_particles(std::vector<Particle> particles) {
    for (const auto& p : particles) {
        if (!finite(p.position) || !finite(p.velocity) ||
            !std::isfinite(p.mass) || p.mass <= 0.0) {
            throw std::invalid_argument("Particle state must be finite and mass positive");
        }
    }
    particles_ = std::move(particles);
}

double LiquidSolver::poly6(double r, double h) {
    if (!(h > 0.0) || r < 0.0 || r >= h) return 0.0;
    const double h2 = h * h;
    const double d = h2 - r * r;
    return (315.0 / (64.0 * pi * std::pow(h, 9))) * d * d * d;
}

Vec3 LiquidSolver::spiky_gradient(Vec3 displacement, double h) {
    const double r = length(displacement);
    if (!(h > 0.0) || r <= 1.0e-12 || r >= h) return {};
    const double scale = -45.0 / (pi * std::pow(h, 6)) * (h - r) * (h - r) / r;
    return displacement * scale;
}

double LiquidSolver::viscosity_laplacian(double r, double h) {
    if (!(h > 0.0) || r < 0.0 || r >= h) return 0.0;
    return 45.0 / (pi * std::pow(h, 6)) * (h - r);
}

void LiquidSolver::compute_density_and_pressure() {
    for (auto& pi_particle : particles_) {
        double density = 0.0;
        for (const auto& pj : particles_) {
            density += pj.mass * poly6(length(pi_particle.position - pj.position),
                                       parameters_.smoothing_length);
        }
        pi_particle.density = std::max(density, density_floor);
        const double raw_pressure = parameters_.sound_speed * parameters_.sound_speed *
                                    (pi_particle.density - parameters_.rest_density);
        pi_particle.pressure = std::max(0.0, raw_pressure);
    }
}

void LiquidSolver::enforce_boundaries(Particle& p) const {
    const double e = parameters_.restitution;
    const auto clamp_axis = [e](double& position, double& velocity, double minimum, double maximum) {
        if (position < minimum) {
            position = minimum;
            if (velocity < 0.0) velocity = -e * velocity;
        } else if (position > maximum) {
            position = maximum;
            if (velocity > 0.0) velocity = -e * velocity;
        }
    };
    clamp_axis(p.position.x, p.velocity.x, parameters_.domain_min.x, parameters_.domain_max.x);
    clamp_axis(p.position.y, p.velocity.y, parameters_.domain_min.y, parameters_.domain_max.y);
    clamp_axis(p.position.z, p.velocity.z, parameters_.domain_min.z, parameters_.domain_max.z);
}

void LiquidSolver::step() {
    compute_density_and_pressure();
    std::vector<Vec3> accelerations(particles_.size(), parameters_.gravity);
    const double h = parameters_.smoothing_length;

    for (std::size_t i = 0; i < particles_.size(); ++i) {
        const auto& pi_particle = particles_[i];
        for (std::size_t j = 0; j < particles_.size(); ++j) {
            if (i == j) continue;
            const auto& pj = particles_[j];
            const Vec3 displacement = pi_particle.position - pj.position;
            const double r = length(displacement);
            if (r >= h || r <= 1.0e-12) continue;

            const Vec3 grad = spiky_gradient(displacement, h);
            const double rhoi = std::max(pi_particle.density, density_floor);
            const double rhoj = std::max(pj.density, density_floor);
            const double pressure_term = pi_particle.pressure / (rhoi * rhoi) +
                                         pj.pressure / (rhoj * rhoj);
            accelerations[i] += grad * (-pj.mass * pressure_term);

            const double visc_term = parameters_.viscosity * pj.mass / rhoj *
                                     viscosity_laplacian(r, h);
            accelerations[i] += (pj.velocity - pi_particle.velocity) * visc_term;
        }
    }

    const double dt = parameters_.time_step;
    for (std::size_t i = 0; i < particles_.size(); ++i) {
        auto& p = particles_[i];
        p.velocity += accelerations[i] * dt;
        p.position += p.velocity * dt;
        enforce_boundaries(p);
        if (!finite(p.position) || !finite(p.velocity) ||
            !std::isfinite(p.density) || !std::isfinite(p.pressure) ||
            !finite(accelerations[i])) {
            throw std::runtime_error("Non-finite state detected in GIFPY 3D liquid solver");
        }
    }
}
} // namespace gifpy
