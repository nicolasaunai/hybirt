#ifndef HYBIRT_MOMENTS_HPP
#define HYBIRT_MOMENTS_HPP

#include "field.hpp"
#include "vecfield.hpp"
#include "population.hpp"

#include <vector>


template<std::size_t dimension>
void density(std::vector<Population<dimension>> const& populations, Field<dimension>& N)
{
    for (auto const& pop : populations)
    {
        for (auto ix = 0; ix < N.data().size(); ++ix)
        {
            N(ix) += pop.density()(ix);
            std::cout << "Density at index " << ix << ": " << N(ix) << "\n";
        }
    }
}

template<std::size_t dimension>
void bulk_velocity(std::vector<Population<dimension>> const& populations, Field<dimension> const& N,
                   VecField<dimension>& V)
{
    for (auto& pop : populations)
    {
        for (auto ix = 0; ix < N.data().size(); ++ix)
        {
            V.x(ix) = pop.flux().x(ix);
            V.y(ix) = pop.flux().y(ix);
            V.z(ix) = pop.flux().z(ix);
        }
    }
    for (auto& pop : populations)
    {
        for (auto ix = 0; ix < N.data().size(); ++ix)
        {
            V.x(ix) /= pop.density()(ix);
            V.y(ix) /= pop.density()(ix);
            V.z(ix) /= pop.density()(ix);
        }
    }
}

#endif
