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
    Faraday(std::shared_ptr<GridLayout<dimension>> grid, double dt)
        : m_grid{grid}
        , m_dt{dt}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& B, VecField<dimension> const& E,
                    VecField<dimension>& Bnew)
    {
        // std::cout << "faraday called\n";
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            for (auto ix = m_grid->ghost_start(Quantity::By, Direction::X);
                 ix <= m_grid->ghost_end(Quantity::By, Direction::X); ++ix)
            {
                // By and Bz loop
                auto const& By = B.y;
                auto const& Bz = B.z;

                auto const& Ey = E.y;
                auto const& Ez = E.z;

                auto& Bnewx = Bnew.x;
                auto& Bnewy = Bnew.y;
                auto& Bnewz = Bnew.z;

                Bnewy(ix) = By(ix) + m_dt * (Ez(ix + 1) - Ez(ix)) / dx;
                Bnewz(ix) = Bz(ix) - m_dt * (Ey(ix + 1) - Ey(ix)) / dx;
                // if (ix == 2)
                // {
                //     std::cout << "Bnewy(" << ix << ") = " << Bnewy(ix) << ", Bnewz(" << ix
                //               << ") = " << Bnewz(ix) << " Ez(ix + 1): " << Ez(ix + 1)
                //               << " Ez(ix): " << Ez(ix) << "\n";
                // }
            }

            for (auto ix = m_grid->ghost_start(Quantity::Bx, Direction::X);
                 ix <= m_grid->ghost_end(Quantity::Bx, Direction::X); ++ix)
            {
                // By and Bz loop
                auto const& Bx = B.x;
                auto& Bnewx    = Bnew.x;
                Bnewx(ix)      = Bx(ix);
            }
        }
        else
            throw std::runtime_error("Faraday not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
    double m_dt;
};

#endif // HYBRIDIR_FARADAY_HPP
