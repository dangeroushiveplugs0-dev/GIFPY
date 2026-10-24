#include "gifpy/liquid_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace gifpy {
namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double density_floor = 1.0e-9;

bool finite(Vec3 v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

struct CellKey {
    std::int64_t x;
    std::int64_t y;
    std::int64_t z;
    bool operator==(const CellKey&) const = default;
};

struct CellKeyHash {
    std::size_t operator()(const CellKey& key) const noexcept {
        const auto mix = [](std::uint64_t value) {
            value ^= value >> 30;
            value *= UINT64_C(0xbf58476d1ce4e5b9);
            value ^= value >> 27;
            value *= UINT64_C(0x94d049bb133111eb);
            return value ^ (value >> 31);
        };
        std::uint64_t hash = mix(static_cast<std::uint64_t>(key.x));
        hash ^= mix(static_cast<std::uint64_t>(key.y) + UINT64_C(0x9e3779b97f4a7c15));
        hash ^= mix(static_cast<std::uint64_t>(key.z) + UINT64_C(0x3c6ef372fe94f82a));
        return static_cast<std::size_t>(hash);
    }
};

std::int64_t cell_coordinate(double coordinate, double cell_width) {
    const long double scaled = std::floor(static_cast<long double>(coordinate) /
                                          static_cast<long double>(cell_width));
    const long double low = static_cast<long double>(std::numeric_limits<std::int64_t>::min()) + 1.0L;
    const long double high = static_cast<long double>(std::numeric_limits<std::int64_t>::max()) - 1.0L;
    if (!std::isfinite(scaled) || scaled < low || scaled > high) {
        throw std::overflow_error("Particle position exceeds spatial-grid coordinate range");
    }
    return static_cast<std::int64_t>(scaled);
}

class SpatialGrid {
public:
    SpatialGrid(const std::vector<Particle>& particles, double cell_width)
        : cell_width_(cell_width) {
        cells_.reserve(particles.size());
        for (std::size_t i = 0; i < particles.size(); ++i) {
            cells_[key_for(particles[i].position)].push_back(i);
        }
    }

    template <typename Function>
    void for_each_candidate(Vec3 position, Function&& function) const {
        const CellKey center = key_for(position);
        for (std::int64_t dx = -1; dx <= 1; ++dx) {
            for (std::int64_t dy = -1; dy <= 1; ++dy) {
                for (std::int64_t dz = -1; dz <= 1; ++dz) {
                    const CellKey key{center.x + dx, center.y + dy, center.z + dz};
                    const auto found = cells_.find(key);
                    if (found == cells_.end()) continue;
                    for (const std::size_t index : found->second) function(index);
                }
            }
        }
    }

private:
    CellKey key_for(Vec3 position) const {
        return {cell_coordinate(position.x, cell_width_),
                cell_coordinate(position.y, cell_width_),
                cell_coordinate(position.z, cell_width_)};
    }

    double cell_width_;
    std::unordered_map<CellKey, std::vector<std::size_t>, CellKeyHash> cells_;
};
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
    const double h = parameters_.smoothing_length;
    const SpatialGrid grid(particles_, h);
    for (auto& pi_particle : particles_) {
        double density = 0.0;
        grid.for_each_candidate(pi_particle.position, [&](std::size_t j) {
            const auto& pj = particles_[j];
            density += pj.mass * poly6(length(pi_particle.position - pj.position), h);
        });
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
    const double h = parameters_.smoothing_length;
    // Particle positions do not change until integration, so one grid can serve
    // both density evaluation and force accumulation during this step.
    const SpatialGrid grid(particles_, h);
    for (auto& pi_particle : particles_) {
        double density = 0.0;
        grid.for_each_candidate(pi_particle.position, [&](std::size_t j) {
            const auto& pj = particles_[j];
            density += pj.mass * poly6(length(pi_particle.position - pj.position), h);
        });
        pi_particle.density = std::max(density, density_floor);
        const double raw_pressure = parameters_.sound_speed * parameters_.sound_speed *
                                    (pi_particle.density - parameters_.rest_density);
        pi_particle.pressure = std::max(0.0, raw_pressure);
    }

    std::vector<Vec3> accelerations(particles_.size(), parameters_.gravity);
    for (std::size_t i = 0; i < particles_.size(); ++i) {
        const auto& pi_particle = particles_[i];
        grid.for_each_candidate(pi_particle.position, [&](std::size_t j) {
            if (i == j) return;
            const auto& pj = particles_[j];
            const Vec3 displacement = pi_particle.position - pj.position;
            const double r = length(displacement);
            if (r >= h || r <= 1.0e-12) return;

            const Vec3 grad = spiky_gradient(displacement, h);
            const double rhoi = std::max(pi_particle.density, density_floor);
            const double rhoj = std::max(pj.density, density_floor);
            const double pressure_term = pi_particle.pressure / (rhoi * rhoi) +
                                         pj.pressure / (rhoj * rhoj);
            accelerations[i] += grad * (-pj.mass * pressure_term);

            const double visc_term = parameters_.viscosity * pj.mass / rhoj *
                                     viscosity_laplacian(r, h);
            accelerations[i] += (pj.velocity - pi_particle.velocity) * visc_term;
        });
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
