#ifndef HYBRIDIR_VECFIELD_HPP
#define HYBRIDIR_VECFIELD_HPP

#include "field.hpp"
#include "gridlayout.hpp"
#include "utils.hpp"

#include <memory>

template<std::size_t dimension>
struct VecField
{
    VecField(std::shared_ptr<GridLayout<dimension>> const& gridlayout,
             std::array<Quantity, 3> quantities)
        : x{gridlayout->allocate(quantities[0])}
        , y{gridlayout->allocate(quantities[1])}
        , z{gridlayout->allocate(quantities[2])}
    {
        // Initialize the vector field with the grid layout
    }
    Field<dimension> x;
    Field<dimension> y;
    Field<dimension> z;
};


#endif // HYBRIDIR_VECFIELD_HPP
