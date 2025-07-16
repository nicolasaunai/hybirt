#ifndef HYBIRT_MOMENTS_HPP
#define HYBIRT_MOMENTS_HPP

#include "field.hpp"
#include "vecfield.hpp"


template<std::size_t dimension>
void bulk_velocity(Field<dimension> const& N, VecField<dimension> const& F, VecField<dimension>& V)
{
    for (auto ix = 0; ix < N.data().size(); ++ix)
    {
        V.x(ix) = F.x(ix) / N(ix);
        V.y(ix) = F.y(ix) / N(ix);
        V.z(ix) = F.z(ix) / N(ix);
    }
}

#endif
