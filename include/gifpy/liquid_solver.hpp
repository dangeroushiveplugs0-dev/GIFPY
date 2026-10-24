#pragma once

#include <vector>

namespace gifpy {

struct Vec3 {
    double x{0.0};
    double y{0.0};
    double z{0.0};
    Vec3& operator+=(const Vec3& other) {
        x += other.x; y += other.y; z += other.z; return *this;
    }
};
Vec3 operator+(Vec3 a, const Vec3& b);
Vec3 operator-(Vec3 a, const Vec3& b);
Vec3 operator*(Vec3 a, double s);
double length(Vec3 v);

struct Particle {
    Vec3 position{};
    Vec3 velocity{};
    double mass{1.0};
    double density{0.0};
    double pressure{0.0};
};

struct Parameters {
    double smoothing_length{0.1};
    double rest_density{1000.0};
    double sound_speed{20.0};
    double viscosity{0.05};
    Vec3 gravity{0.0, -9.81, 0.0};
    double time_step{0.0005};
    Vec3 domain_min{0.0, 0.0, 0.0};
    Vec3 domain_max{1.0, 1.0, 1.0};
    double restitution{0.0};
};

class LiquidSolver {
public:
    explicit LiquidSolver(Parameters parameters = {});
    void set_particles(std::vector<Particle> particles);
    const std::vector<Particle>& particles() const noexcept { return particles_; }
    const Parameters& parameters() const noexcept { return parameters_; }
    void compute_density_and_pressure();
    void step();

    static double poly6(double distance, double h);
    static Vec3 spiky_gradient(Vec3 displacement, double h);
    static double viscosity_laplacian(double distance, double h);

private:
    Parameters parameters_;
    std::vector<Particle> particles_;
    void validate_parameters() const;
    void enforce_boundaries(Particle& particle) const;
};

} // namespace gifpy
