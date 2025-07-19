#include "vecfield.hpp"
#include "field.hpp"

#include "faraday.hpp"
#include "ampere.hpp"
#include "ohm.hpp"
#include "utils.hpp"
#include "gridlayout.hpp"
#include "boundary_condition.hpp"
#include "moments.hpp"
#include "pusher.hpp"
#include "diagnostics.hpp"
#include "population.hpp"

#include "highfive/highfive.hpp"

#include <iostream>
#include <array>
#include <vector>
#include <cstdint>
#include <memory>
#include <algorithm>




template<std::size_t dimension>
void average(Field<dimension> const& F1, Field<dimension> const& F2, Field<dimension>& Favg)
{
    std::transform(F1.begin(), F1.end(), F2.begin(), Favg.begin(),
                   [](double const& f1, double const& f2) { return 0.5 * (f1 + f2); });
}


template<std::size_t dimension>
void average(VecField<dimension> const& V1, VecField<dimension> const& V2,
             VecField<dimension>& Vavg)
{
    average(V1.x, V2.x, Vavg.x);
    average(V1.y, V2.y, Vavg.y);
    average(V1.z, V2.z, Vavg.z);
}


double bx(double x)
{
    // Placeholder for a function that returns Bx based on x
    return 0.0; // Example value
}

double by(double x)
{
    // Placeholder for a function that returns By based on x
    return 1.0; // Example value
}


double bz(double x)
{
    // Placeholder for a function that returns Bz based on x
    return 0.0; // Example value
}

double density(double x)
{
    // Placeholder for a function that returns density based on x
    return 1.0; // Example value
}


void magnetic_init(VecField<1>& B, GridLayout<1> const& layout)
{
    // Initialize magnetic field B
    for (auto ix = layout.primal_dom_start(Direction::X); ix <= layout.primal_dom_end(Direction::X);
         ++ix)
    {
        auto x = layout.coordinate(Direction::X, Quantity::Bx, ix);

        B.x(ix) = bx(x); // Bx
    }
    for (auto ix = layout.dual_dom_start(Direction::X); ix <= layout.dual_dom_end(Direction::X);
         ++ix)
    {
        auto x = layout.coordinate(Direction::X, Quantity::By, ix);

        B.y(ix) = by(x); // By
        B.z(ix) = bz(x); // Bz, uniform magnetic field in z-direction
    }
}




int main()
{
    double time                     = 0.;
    double final_time               = 4.0000;
    double dt                       = 0.0001;
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {100};
    std::array<double, dimension> cell_size      = {0.2};
    auto constexpr nbr_ghosts                    = 1;
    auto constexpr nppc                          = 100;

    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> E{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> Enew{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> Bnew{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> Eavg{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> Bavg{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};
    VecField<dimension> J{layout, {Quantity::Jx, Quantity::Jy, Quantity::Jz}};
    VecField<dimension> V{layout, {Quantity::Vx, Quantity::Vy, Quantity::Vz}};
    Field<dimension> N{layout->allocate(Quantity::N), Quantity::N};

    auto boundary_condition = BoundaryConditionFactory<dimension>::create("periodic", layout);

    std::vector<Population<1>> populations;
    populations.emplace_back("main", layout);
    for (auto& pop : populations)
        pop.load_particles(nppc, density);


    magnetic_init(B, *layout);
    boundary_condition->fill(B);

    Faraday<dimension> faraday{layout, dt};
    Ampere<dimension> ampere{layout};
    Ohm<dimension> ohm{layout};
    Boris<dimension> push{layout, dt};



    ampere(B, J);
    boundary_condition->fill(J);
    for (auto& pop : populations)
    {
        pop.deposit();
        boundary_condition->fill(pop.flux());
        boundary_condition->fill(pop.density());
    }

    total_density(populations, N);
    bulk_velocity<dimension>(populations, N, V);
    ohm(B, J, N, V, E);
    boundary_condition->fill(E);

    diags_write_fields(B, E, V, N, time, HighFive::File::Truncate);
    diags_write_particles(populations, time, HighFive::File::Truncate);

    while (time < final_time)
    {
        std::cout << "Time: " << time << " / " << final_time << "\n";



        // predictor 1
        std::cout << "Predictor 1 step-------------\n";
        faraday(B, E, Bnew);
        ampere(Bnew, J);
        boundary_condition->fill(J);
        ohm(Bnew, J, N, V, Enew);
        average(E, Enew, Eavg);
        average(B, Bnew, Bavg);
        boundary_condition->fill(Eavg);
        boundary_condition->fill(Bavg);




        // save particles at t=n before updating to predicted t=n+1
        // because predictor 2 will re-start from t=n
        auto save{populations};
        for (auto& pop : save)
        {
            push(pop.particles(), Eavg, Bavg);
            boundary_condition->particles(pop.particles());
            pop.deposit();
            boundary_condition->fill(pop.flux());
            boundary_condition->fill(pop.density());
        }
        total_density(populations, N);
        bulk_velocity<dimension>(populations, N, V);
        save.clear(); // clear the save vector to free memory




        // predictor 2
        std::cout << "Predictor 2 step-------------\n";
        faraday(B, Eavg, Bnew);
        ampere(Bnew, J);
        boundary_condition->fill(J);
        ohm(Bnew, J, N, V, Enew);
        average(E, Enew, Eavg);
        average(B, Bnew, Bavg);
        boundary_condition->fill(Eavg);
        boundary_condition->fill(Bavg);




        // populations = save;
        for (auto& pop : populations)
        {
            push(pop.particles(), Eavg, Bavg);
            boundary_condition->particles(pop.particles());
            pop.deposit();
            boundary_condition->fill(pop.flux());
            boundary_condition->fill(pop.density());
        }
        total_density(populations, N);
        bulk_velocity<dimension>(populations, N, V);




        // corrector
        std::cout << "Corrector step-------------\n";
        faraday(B, Eavg, B);
        ampere(B, J);
        boundary_condition->fill(J);
        ohm(B, J, N, V, E);
        boundary_condition->fill(E);


        time += dt;
        diags_write_fields(B, E, V, N, time);
        std::cout << "**********************************\n";
        // diags_write_particles(populations, time);
    }


    return 0;
}
