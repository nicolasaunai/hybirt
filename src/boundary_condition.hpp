#ifndef HYBIRT_BOUNDARY_CONDITION_HPP
#define HYBIRT_BOUNDARY_CONDITION_HPP

#include "field.hpp"
#include "vecfield.hpp"
#include "pusher.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <stdexcept>


template<std::size_t dimension>
class BoundaryCondition
{
public:
    BoundaryCondition(std::shared_ptr<GridLayout<dimension>> const& grid)
        : m_grid(grid)
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    virtual void fill(Field<dimension>& field) = 0;

    void fill(VecField<dimension>& vecfield)
    {
        fill(vecfield.x);
        fill(vecfield.y);
        fill(vecfield.z);
    }

    virtual void particles(std::vector<Particle<dimension>>& particles) = 0;

protected:
    std::shared_ptr<GridLayout<dimension>> const& m_grid;
};


template<std::size_t dimension>
class PeriodicBoundaryCondition : public BoundaryCondition<dimension>
{
public:
    PeriodicBoundaryCondition(std::shared_ptr<GridLayout<dimension>> const& grid)
        : BoundaryCondition<dimension>(grid)
    {
    }

    void fill(Field<dimension>& field) override
    {
        if constexpr (dimension == 1)
        {
            // left side
            auto gsi = this->m_grid->ghost_start(field.quantity(), Direction::X);
            auto dsi = this->m_grid->dom_start(field.quantity(), Direction::X);
            auto dei = this->m_grid->dom_end(field.quantity(), Direction::X);
            auto gei = this->m_grid->ghost_end(field.quantity(), Direction::X);

            auto const dom_size = this->m_grid->nbr_cells(Direction::X);


            // std::cout << "Filling left side\n";
            for (auto ix_left = gsi; ix_left < dsi; ++ix_left)
            {
                auto const ix_right = ix_left + dom_size;
                // std::cout << "field(" << ix_left << ") = field(" << ix_right << ")\n";
                field(ix_left) = field(ix_right);
            }

            // std::cout << "Filling right side\n";
            for (auto ix_right = gei; ix_right > dei; --ix_right)
            {
                auto const ix_left = ix_right - dom_size;
                // std::cout << "field(" << ix_right << ") = field(" << ix_left << ")\n";
                field(ix_right) = field(ix_left);
            }
        }
    }

    void particles(std::vector<Particle<dimension>>& particles) override
    {
        if constexpr (dimension == 1)
        {
            for (auto& particle : particles)
            {
                double cell = static_cast<int>(particle.position[0]
                                               / this->m_grid->cell_size(Direction::X));

                // particles left the right border injected on left side
                if (cell > this->m_grid->dual_dom_end(Direction::X))
                {
                    cell -= this->m_grid->nbr_cells(Direction::X);
                    particle.position[0] = cell;
                }
                // particles left the left border injected on right side
                else if (cell < this->m_grid->dual_dom_start(Direction::X))
                {
                    // Wrap around tothe right side
                    cell += this->m_grid->nbr_cells(Direction::X);
                    particle.position[0] = cell;
                }
            }
        }
    }
};



template<std::size_t dimension>
class BoundaryConditionFactory
{
public:
    static std::unique_ptr<BoundaryCondition<dimension>>
    create(std::string const& type, std::shared_ptr<GridLayout<dimension>> grid)
    {
        if (type == "periodic")
            return std::make_unique<PeriodicBoundaryCondition<dimension>>(grid);
        // Add more boundary condition types as needed
        throw std::runtime_error("Unknown boundary condition type: " + type);
    }
};



#endif
