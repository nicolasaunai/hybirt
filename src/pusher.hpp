#ifndef PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP
#define PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP


#include "vecfield.hpp"
#include "particle.hpp"

#include <cstddef>
#include <vector>


template<std::size_t dimension>
class Pusher
{
protected:
    std::shared_ptr<GridLayout<dimension>> layout_;
    double dt_;

public:
    Pusher(std::shared_ptr<GridLayout<dimension>> layout, double dt)
        : layout_(layout)
        , dt_(dt)
    {
    }

    virtual void operator()(std::vector<Particle<dimension>>& particles,
                            VecField<dimension> const& E, VecField<dimension> const& B)
        = 0;

    virtual ~Pusher() {}
};



template<std::size_t dimension>
class Boris : public Pusher<dimension>
{
public:
    Boris(std::shared_ptr<GridLayout<dimension>> layout, double dt)
        : Pusher<dimension>{layout, dt}
    {
    }

    void operator()(std::vector<Particle<dimension>>& particles, VecField<dimension> const& E,
                    VecField<dimension> const& B) override
    {
        for (auto& particle : particles)
        {
            // half step push position from t=n to t=n+1/2
            // needed because fields are defined at t=n+1/2
            for (auto dim = 0; dim < dimension; ++dim)
            {
                particle.position[dim] += particle.v[dim] * this->dt_ * 0.5;
            }

            double const iCell_float = particle.position[0] / this->layout_->cell_size(Direction::X)
                                       + this->layout_->dual_dom_start(Direction::X);
            int const iCell       = static_cast<int>(iCell_float);
            double const reminder = iCell_float - iCell;
            double const qdto2m   = particle.charge * this->dt_ / (2.0 * particle.mass);


            // Interpolate E and B fields at particle position
            auto const dx   = this->layout_->cell_size(Direction::X);
            double const ex = interpolate(E.x, iCell, reminder);
            double const ey = interpolate(E.y, iCell, reminder);
            double const ez = interpolate(E.z, iCell, reminder);
            double const bx = interpolate(B.x, iCell, reminder);
            double const by = interpolate(B.y, iCell, reminder);
            double const bz = interpolate(B.z, iCell, reminder);


            // Calculate the half-step velocity
            auto const vminus_x = particle.v[0] + qdto2m * ex;
            auto const vminus_y = particle.v[1] + qdto2m * ey;
            auto const vminus_z = particle.v[2] + qdto2m * ez;



            auto const vprime_x = vminus_x + qdto2m * (vminus_y * bz - vminus_z * by);
            auto const vprime_y = vminus_y + qdto2m * (vminus_z * bx - vminus_x * bz);
            auto const vprime_z = vminus_z + qdto2m * (vminus_x * by - vminus_y * bx);


            auto const s       = 2 * qdto2m / (1 + qdto2m * qdto2m * (bx * bx + by * by + bz * bz));
            auto const vplus_x = vminus_x + s * (vprime_y * bz - vprime_z * by);
            auto const vplus_y = vminus_y + s * (vprime_z * bx - vprime_x * bz);
            auto const vplus_z = vminus_z + s * (vprime_x * by - vprime_y * bx);

            particle.v[0] = vplus_x + qdto2m * ex;
            particle.v[1] = vplus_y + qdto2m * ey;
            particle.v[2] = vplus_z + qdto2m * ez;

            // velocity is now at t=n+1
            // so we can compute the other half of the position update
            for (auto dim = 0; dim < dimension; ++dim)
            {
                particle.position[dim] += particle.v[dim] * this->dt_ * 0.5;
            }
        }
    }

private:
    double interpolate(Field<dimension> const& field, int iCell, double reminder) const
    {
        return field(iCell) * (1.0 - reminder) + field(iCell + 1) * reminder;
    }
};


#endif
