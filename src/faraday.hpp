#ifndef HYBRIDIR_FARADAY_HPP
#define HYBRIDIR_FARADAY_HPP

#include "vecfield.hpp"
#include "gridlayout.hpp"
#include "utils.hpp"

#include <cstddef>
#include <iostream>
#include <memory>

template<std::size_t dimension>
class Faraday
{
public:
    void operator()(VecField<dimension> const& B, VecField<dimension> const& E,
                    VecField<dimension>& Bnew)
    {
        std::cout << "faraday called\n";
        for (auto ix = m_grid->dual_dom_start(Direction::X);
             ix <= m_grid->dual_dom_end(Direction::X); ++ix)
        {
            // By and Bz loop
            auto const& By = B.y;
            auto const& Bz = B.z;

            auto const& Ey = E.y;
            auto const& Ez = E.z;

            auto const& Bnewy = Bnew.y;
            auto const& Bnewz = Bnew.z;
        }
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};

#endif // HYBRIDIR_FARADAY_HPP
