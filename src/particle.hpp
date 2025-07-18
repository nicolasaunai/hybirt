#ifndef HYBIRT_PARTICLE_HPP
#define HYBIRT_PARTICLE_HPP

#include <array>
#include <cstddef>

template<std::size_t dimension>
struct Particle
{
    std::array<double, dimension> position;
    std::array<double, 3> v; // velocity
    double weight;
    double mass;
    double charge;
};
#endif // HYBIRT_PARTICLE_HPP
