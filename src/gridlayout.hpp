#ifndef HYBIRT_GRIDLAYOUT_HPP
#define HYBIRT_GRIDLAYOUT_HPP

#include "utils.hpp"

#include <cstddef>
#include <array>
#include <stdexcept>


enum class Quantity { E, Ex, Ey, Ez, B, Bx, By, Bz, J, Jx, Jy, Jz, N, V, Vx, Vy, Vz };


template<std::size_t dimension>
class GridLayout
{
public:
    GridLayout(std::array<std::size_t, dimension> nbr_cells,
               std::array<double, dimension> cell_size, std::size_t nbr_ghosts)
        : m_nbr_cells{nbr_cells}
        , m_cell_size{cell_size}
        , m_nbr_ghosts{nbr_ghosts}
    {
    }


    auto dual_dom_start(Direction dir_idx) const { return m_nbr_ghosts; }
    auto dual_dom_end(Direction dir_idx) const
    {
        return m_nbr_cells[dir_idx] + dual_dom_start(dir_idx) - 1;
    }

    auto primal_dom_start(Direction dir_idx) const { return m_nbr_ghosts; }
    auto primal_dom_end(Direction dir_idx) const
    {
        return m_nbr_cells[dir_idx] + primal_dom_start(dir_idx);
    }

    auto cell_size(Direction dir_idx) const { return m_cell_size[dir_idx]; }

    auto allocate(Quantity qty)
    {
        // Placeholder for allocation logic
        // This would typically allocate memory for the field based on the quantity and direction
        auto centering = centerings(qty);
        if constexpr (dimension == 1)
            return std::array<std::size_t, dimension>{m_nbr_cells[0] + 2 * m_nbr_ghosts
                                                      + centering[0]};
        else if constexpr (dimension == 2)
            return std::array<std::size_t, dimension>{
                m_nbr_cells[0] + 2 * m_nbr_ghosts + centering[0],
                m_nbr_cells[1] + 2 * m_nbr_ghosts + centering[1]};

        else if constexpr (dimension == 3)
            return std::array<std::size_t, dimension>{
                m_nbr_cells[0] + 2 * m_nbr_ghosts + centering[0],
                m_nbr_cells[1] + 2 * m_nbr_ghosts + centering[1],
                m_nbr_cells[2] + 2 * m_nbr_ghosts + centering[2]};
        else
            throw std::runtime_error("Unsupported dimension");
    }

private:
    std::array<std::size_t, dimension> centerings(Quantity qty)
    {
        switch (qty)
        {
            case Quantity::Jx:
            case Quantity::Ex:
                if constexpr (dimension == 1)
                    return {0};
                else if constexpr (dimension == 2)
                    return {0, 1};
                else if constexpr (dimension == 3)
                    return {0, 1, 1};

            case Quantity::Jy:
            case Quantity::Ey:
                if constexpr (dimension == 1)
                    return {1};
                else if constexpr (dimension == 2)
                    return {1, 0};
                else if constexpr (dimension == 3)
                    return {1, 0, 1};

            case Quantity::Jz:
            case Quantity::Ez:
                if constexpr (dimension == 1)
                    return {1};
                else if constexpr (dimension == 2)
                    return {1, 1};
                else if constexpr (dimension == 3)
                    return {1, 1, 0};

            case Quantity::Bx:
                if constexpr (dimension == 1)
                    return {1};
                else if constexpr (dimension == 2)
                    return {1, 0};
                else if constexpr (dimension == 3)
                    return {1, 0, 0};

            case Quantity::By:
                if constexpr (dimension == 1)
                    return {0};
                else if constexpr (dimension == 2)
                    return {0, 1};
                else if constexpr (dimension == 3)
                    return {0, 1, 0};

            case Quantity::Bz:
                if constexpr (dimension == 1)
                    return {0};
                else if constexpr (dimension == 2)
                    return {0, 0};
                else if constexpr (dimension == 3)
                    return {0, 0, 1};

            case Quantity::N:
            case Quantity::V:
                if constexpr (dimension == 1)
                    return {1};
                else if constexpr (dimension == 2)
                    return {1, 1};
                else if constexpr (dimension == 3)
                    return {1, 1, 1};

            default: throw std::runtime_error{"Unknown quantity"};
        }
    }


    std::array<std::size_t, dimension> m_nbr_cells;
    std::array<double, dimension> m_cell_size;
    std::size_t m_nbr_ghosts;
};

#endif // HYBIRT_GRIDLAYOUT_HPP
