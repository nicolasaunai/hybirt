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
private:
    static constexpr auto dual   = 0;
    static constexpr auto primal = 1;

public:
    GridLayout(std::array<std::size_t, dimension> nbr_cells,
               std::array<double, dimension> cell_size, std::size_t nbr_ghosts)
        : m_nbr_cells{nbr_cells}
        , m_cell_size{cell_size}
        , m_nbr_ghosts{nbr_ghosts}
    {
    }

    auto nbr_cells(Direction dir_idx) const { return m_nbr_cells[dir_idx]; }

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

    auto ghost_start(Quantity qty, Direction dir_idx) const { return 0; }

    auto ghost_end(Quantity qty, Direction dir_idx) const { return allocate(qty)[dir_idx] - 1; }

    auto dom_start(Quantity qty, Direction dir_idx) const
    {
        auto const centering = centerings(qty);
        if constexpr (dimension == 1)
        {
            return centering[0] == dual ? dual_dom_start(dir_idx) : primal_dom_start(dir_idx);
        }
    }
    auto dom_end(Quantity qty, Direction dir_idx) const
    {
        auto const centering = centerings(qty);
        if constexpr (dimension == 1)
        {
            return centering[0] == dual ? dual_dom_end(dir_idx) : primal_dom_end(dir_idx);
        }
    }

    auto cell_size(Direction dir_idx) const { return m_cell_size[dir_idx]; }

    template<typename... Index>
    auto coordinate(Direction dir, Quantity qty, Index... index) const
    {
        return 0.; // Placeholder for coordinate calculation logic
    }

    auto allocate(Quantity qty) const
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
    std::array<std::size_t, dimension> centerings(Quantity qty) const
    {
        switch (qty)
        {
            case Quantity::Jx:
            case Quantity::Ex:
                if constexpr (dimension == 1)
                    return {dual};
                else if constexpr (dimension == 2)
                    return {dual, primal};
                else if constexpr (dimension == 3)
                    return {dual, primal, primal};

            case Quantity::Jy:
            case Quantity::Ey:
                if constexpr (dimension == 1)
                    return {primal};
                else if constexpr (dimension == 2)
                    return {primal, dual};
                else if constexpr (dimension == 3)
                    return {primal, dual, primal};

            case Quantity::Jz:
            case Quantity::Ez:
                if constexpr (dimension == 1)
                    return {primal};
                else if constexpr (dimension == 2)
                    return {primal, primal};
                else if constexpr (dimension == 3)
                    return {primal, primal, dual};

            case Quantity::Bx:
                if constexpr (dimension == 1)
                    return {primal};
                else if constexpr (dimension == 2)
                    return {primal, dual};
                else if constexpr (dimension == 3)
                    return {primal, dual, dual};

            case Quantity::By:
                if constexpr (dimension == 1)
                    return {dual};
                else if constexpr (dimension == 2)
                    return {dual, primal};
                else if constexpr (dimension == 3)
                    return {dual, primal, dual};

            case Quantity::Bz:
                if constexpr (dimension == 1)
                    return {dual};
                else if constexpr (dimension == 2)
                    return {dual, dual};
                else if constexpr (dimension == 3)
                    return {dual, dual, primal};

            case Quantity::N:
            case Quantity::V:
            case Quantity::Vx:
            case Quantity::Vy:
            case Quantity::Vz:
                if constexpr (dimension == 1)
                    return {primal};
                else if constexpr (dimension == 2)
                    return {primal};
                else if constexpr (dimension == 3)
                    return {primal, primal, primal};

            default: throw std::runtime_error{"Unknown quantity"};
        }
    }


    std::array<std::size_t, dimension> m_nbr_cells;
    std::array<double, dimension> m_cell_size;
    std::size_t m_nbr_ghosts;
};

#endif // HYBIRT_GRIDLAYOUT_HPP
