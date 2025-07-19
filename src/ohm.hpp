#ifndef HYBRIDIR_OHM_HPP
#define HYBRIDIR_OHM_HPP

#include "vecfield.hpp"

#include <cstddef>
#include <memory>
#include <iostream>

template<std::size_t dimension>
class Ohm
{
public:
    Ohm(std::shared_ptr<GridLayout<dimension>> grid)
        : m_grid{grid}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }
    void operator()(VecField<dimension> const& B, VecField<dimension> const& J, Field<dimension>& N,
                    VecField<dimension> const& V, VecField<dimension>& Enew)

    {
        // Placeholder implementation
        // std::cout << "ohm called\n";
        // Ex = -(Vy * Bz - Vz * By) + (JyBz - JzBy) / N;

        auto const dx = m_grid->cell_size(Direction::X);
        if constexpr (dimension == 1)
        {
            // Ex is dual in x
            for (auto ix = m_grid->dual_dom_start(Direction::X);
                 ix <= m_grid->dual_dom_end(Direction::X); ++ix)
            {
                auto const& By = B.y;
                auto const& Bz = B.z;

                auto const& Jx = J.x;
                auto const& Jy = J.y;
                auto const& Jz = J.z;

                auto const& Vy = V.y;
                auto const& Vz = V.z;

                auto& Ex = Enew.x;

                auto const Vy_dual = 0.5 * (Vy(ix + 1) + Vy(ix));
                auto const Vz_dual = 0.5 * (Vz(ix + 1) + Vz(ix));
                auto const N_dual  = 0.5 * (N(ix + 1) + N(ix));
                auto const Jy_dual = 0.5 * (Jy(ix + 1) + Jy(ix));
                auto const Jz_dual = 0.5 * (Jz(ix + 1) + Jz(ix));

                auto const ideal_x = -(Vy_dual * Bz(ix) - Vz_dual * By(ix));
                auto const hall_x  = (Jy_dual * Bz(ix) - Jz_dual * By(ix)) / N_dual;


                // if (ix == 2)
                //     std::cout << "hall_x at ix=" << ix << ": " << hall_x << " jy: " << Jy(ix)
                //               << ", jz: " << Jz(ix) << ", by: " << By(ix) << ", bz: " << Bz(ix)
                //               << ", n: " << N_dual << " N(ix + 1): " << N(ix + 1)
                //               << " N(ix) = " << N(ix) << "\n";

                Ex(ix) = ideal_x + 1 * hall_x + 0.000 * Jx(ix);
                // - 0.001 * (Jx(ix + 1) - 2 * Jx(ix) + Jx(ix - 1))
                //       / (dx * dx); // + 0.001*Jx(ix);
            }

            // Ey is primal in x, so is Ez
            // Ey = -(Vz * Bx - Vx * Bz) + (JzBx - JxBz) / N;
            // Ez = -(Vx * By - Vy * Bx) + (JxBy - JyBx) / N;
            for (auto ix = m_grid->primal_dom_start(Direction::X);
                 ix <= m_grid->primal_dom_end(Direction::X); ++ix)
            {
                auto const& Bx = B.x;
                auto const& By = B.y;
                auto const& Bz = B.z;

                auto const& Jx = J.x;
                auto const& Jy = J.y;
                auto const& Jz = J.z;

                auto const& Vx = V.x;
                auto const& Vy = V.y;
                auto const& Vz = V.z;

                auto& Ey = Enew.y;
                auto& Ez = Enew.z;

                auto const Jx_primal = 0.5 * (Jx(ix) + Jx(ix - 1));
                auto const Bz_primal = 0.5 * (Bz(ix) + Bz(ix - 1));
                auto const By_primal = 0.5 * (By(ix) + By(ix - 1));

                auto const ideal_y = -(Vz(ix) * Bx(ix) - Vx(ix) * Bz_primal);
                auto const hall_y  = (Jz(ix) * Bx(ix) - Jx_primal * Bz_primal) / N(ix);

                auto const ideal_z = -(Vx(ix) * By_primal - Vy(ix) * Bx(ix));
                auto const hall_z  = (Jx_primal * By_primal - Jy(ix) * Bx(ix)) / N(ix);

                Ey(ix) = ideal_y + 1 * hall_y + 0.000 * Jy(ix);
                Ez(ix) = ideal_z + 1 * hall_z + 0.000 * Jz(ix);
            }
        }
        else
        {
            throw std::runtime_error("Ohm's law not implemented for this dimension");
        }
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};

#endif // HYBRIDIR_OHM_HPP
