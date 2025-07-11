#ifndef HYBRIDIR_OHM_HPP
#define HYBRIDIR_OHM_HPP

#include "vecfield.hpp"

#include <cstddef>
#include <iostream>

template<std::size_t dimension>
class Ohm
{
public:
    void operator()(VecField<dimension> const& B, VecField<dimension> const& J, Field<dimension>& N,
                    VecField<dimension>& V, VecField<dimension>& Enew)

    {
        // Placeholder implementation
        std::cout << "ohm called\n";
    }
};

#endif // HYBRIDIR_OHM_HPP
