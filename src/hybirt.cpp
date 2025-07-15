#include "vecfield.hpp"
#include "field.hpp"

#include "faraday.hpp"
#include "ampere.hpp"
#include "ohm.hpp"
#include "utils.hpp"
#include "gridlayout.hpp"
#include "boundary_condition.hpp"
#include "moments.hpp"

#include <iostream>
#include <array>
#include <vector>
#include <cstdint>
#include <memory>




template<std::size_t dimension>
struct Particle
{
    std::array<double, dimension> position;
    std::array<double, 3> v; // velocity
};




template<std::size_t dimension>
void push(std::vector<Particle<1>>& particles, VecField<dimension> const& Eavg,
          VecField<dimension> const& Bavg, double dt)
{
    // Placeholder implementation for pushing particles
}


template<std::size_t dimension>
void deposit(std::vector<Particle<1>>& particles, Field<dimension>& N, VecField<dimension>& V)
{
}



template<std::size_t dimension>
void average(VecField<dimension> const& V1, VecField<dimension> const& V2,
             VecField<dimension>& Vavg)
{
    // placeholder
}


int main()
{
    double time                     = 0.;
    double final_time               = 0.004;
    double dt                       = 0.001;
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {10};
    std::array<double, dimension> cell_size      = {0.1};
    auto constexpr nbr_ghosts                    = 1;
    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> E{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> Enew{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> Bnew{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> Eavg{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> Bavg{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> J{layout, {Quantity::Jx, Quantity::Jy, Quantity::Jz}};
    VecField<dimension> F{layout, {Quantity::Vx, Quantity::Vy, Quantity::Vz}};
    VecField<dimension> V{layout, {Quantity::Vx, Quantity::Vy, Quantity::Vz}};
    Field<dimension> N{layout->allocate(Quantity::N), Quantity::N};
    std::vector<Particle<dimension>> particles;


    Faraday<dimension> faraday{layout, dt};
    Ampere<dimension> ampere{layout};
    Ohm<dimension> ohm;

    auto boundary_condition = BoundaryConditionFactory<dimension>::create("periodic", layout);

    while (time < final_time)
    {
        // predictor 1
        faraday(B, E, Bnew);
        ampere(B, J);
        boundary_condition->fill(J);
        ohm(B, J, N, V, Enew);

        average(E, Enew, Eavg);
        average(B, Bnew, Bavg);
        boundary_condition->fill(Eavg);
        boundary_condition->fill(Bavg);

        push(particles, Eavg, Bavg, dt);
        deposit(particles, N, F);
        bulk_velocity<dimension>(N, F, V);
        boundary_condition->fill(N);
        boundary_condition->fill(V);

        // predictor 2
        faraday(B, Eavg, Bnew);
        ampere(Bnew, J);
        boundary_condition->fill(J);
        ohm(Bnew, J, N, V, Enew);

        average(E, Enew, Eavg);
        average(B, Bnew, Bavg);
        boundary_condition->fill(Eavg);
        boundary_condition->fill(Bavg);

        push(particles, Eavg, Bavg, dt);
        deposit(particles, N, V);
        deposit(particles, N, F);
        bulk_velocity<dimension>(N, F, V);
        boundary_condition->fill(N);
        boundary_condition->fill(V);

        // corrector
        faraday(B, Eavg, B);
        ampere(B, J);
        boundary_condition->fill(J);
        ohm(B, J, N, V, E);
        boundary_condition->fill(E);


        time += dt;
    }


    return 0;
}
