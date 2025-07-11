#ifndef HYBRIDIR_AMPERE_HPP
#define HYBRIDIR_AMPERE_HPP

#include "vecfield.hpp"

#include <cstddef>
#include <iostream>

template<std::size_t dimension>
class Ampere
{
public:
    void operator()(VecField<dimension> const& B, VecField<dimension>& J)
    {
        // Placeholder implementation
        std::cout << "ampere called\n";
    }
};

#endif // HYBRIDIR_AMPERE_HPP
