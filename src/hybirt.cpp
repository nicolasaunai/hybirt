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

#include "highfive/highfive.hpp"


#include <iostream>
#include <array>
#include <vector>
#include <cstdint>
#include <memory>
#include <algorithm>
#include <random>
#include <optional>
#include <iomanip>




template<std::size_t dimension>
void deposit(std::vector<Particle<1>>& particles, Field<dimension>& N, VecField<dimension>& F,
             GridLayout<dimension> const& grid)
{
    for (auto& n : N.data())
    {
        n = 0.0; // Reset the field
    }

    for (auto& fx : F.x.data())
        fx = 0.0;
    for (auto& fy : F.y.data())
        fy = 0.0;
    for (auto& fz : F.z.data())
        fz = 0.0;

    for (auto const& particle : particles)
    {
        double const iCell_float = particle.position[0] / grid.cell_size(Direction::X);
        int const iCell_         = static_cast<int>(iCell_float);
        double const reminder    = iCell_float - iCell_;
        auto const iCell         = iCell_ + grid.dual_dom_start(Direction::X);


        std::cout << "Depositing particle position: " << particle.position[0]
                  << ", iCell: " << iCell << ", reminder: " << reminder
                  << ", weight: " << particle.weight << "\n";

        N(iCell) += particle.weight * (1.0 - reminder);
        N(iCell + 1) += particle.weight * reminder;

        F.x(iCell) += particle.weight * (1.0 - reminder) * particle.v[0];
        F.x(iCell + 1) += particle.weight * reminder * particle.v[0];

        F.y(iCell) += particle.weight * (1.0 - reminder) * particle.v[1];
        F.y(iCell + 1) += particle.weight * reminder * particle.v[1];

        F.z(iCell) += particle.weight * (1.0 - reminder) * particle.v[2];
        F.z(iCell + 1) += particle.weight * reminder * particle.v[2];
    }
}


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
    return 1.0; // Example value
}

double by(double x)
{
    // Placeholder for a function that returns By based on x
    return 0.0; // Example value
}


double bz(double x)
{
    // Placeholder for a function that returns Bz based on x
    return 1.0; // Example value
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
        B.y(ix) = by(x); // By
        B.z(ix) = bz(x); // Bz, uniform magnetic field in z-direction
    }
}

std::mt19937_64 getRNG(std::optional<std::size_t> const& seed)
{
    if (!seed.has_value())
    {
        std::random_device randSeed;
        std::seed_seq seed_seq{randSeed(), randSeed(), randSeed(), randSeed(),
                               randSeed(), randSeed(), randSeed(), randSeed()};
        return std::mt19937_64(seed_seq);
    }
    return std::mt19937_64(*seed);
}


void maxwellianVelocity(std::array<double, 3> const& V, std::array<double, 3> const& Vth,
                        std::mt19937_64& generator, std::array<double, 3>& partVelocity)
{
    std::normal_distribution<> maxwellX(V[0], Vth[0]);
    std::normal_distribution<> maxwellY(V[1], Vth[1]);
    std::normal_distribution<> maxwellZ(V[2], Vth[2]);

    partVelocity[0] = maxwellX(generator);
    partVelocity[1] = maxwellY(generator);
    partVelocity[2] = maxwellZ(generator);
}


template<typename T>
auto to_string_fixed_width(T const& value, std::size_t const& precision, std::size_t const& width,
                           char const& fill = '0')
{
    std::ostringstream out;
    out.width(width);
    out.precision(precision);
    out << std::setfill(fill) << std::fixed << value;
    return out.str();
}

void load_particles(std::vector<Particle<1>>& particles, GridLayout<1> const& layout, int nppc)
{
    auto randGen = getRNG(std::nullopt);
    std::array<double, 3> Vth{0.1, 0.1, 0.1}; // thermal velocity in each direction
    std::array<double, 3> V{1.0, 0.0, 0.0};   // bulk velocity

    for (auto iCell = layout.dual_dom_start(Direction::X);
         iCell <= layout.dual_dom_end(Direction::X); ++iCell)
    {
        auto const x      = layout.cell_coordinate(Direction::X, iCell);
        auto cell_density = density(x);

        auto cell_weight = cell_density / nppc;
        for (auto partIdx = 0; partIdx < nppc; ++partIdx)
        {
            Particle<1> particle;
            particle.position[0] = x + 0.0 * layout.cell_size(Direction::X); // center of the cell
            maxwellianVelocity(V, Vth, randGen, particle.v);
            particle.weight = cell_weight;
            particle.mass   = 1.0; // hard-coded mass
            particle.charge = 1.0; // hard-coded charge

            std::cout << "Particle " << particles.size() << ": "
                      << "Position: " << particle.position[0] << ", "
                      << "Velocity: (" << particle.v[0] << ", " << particle.v[1] << ", "
                      << particle.v[2] << "), "
                      << "Weight: " << particle.weight << "\n";
            particles.push_back(particle);
        }
    }
}




int main()
{
    double time                     = 0.;
    double final_time               = 0.001;
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

    auto constexpr nppc = 1;
    std::vector<Particle<dimension>> particles;


    magnetic_init(B, *layout);
    load_particles(particles, *layout, nppc);
    std::cout << "Loaded " << particles.size() << " particles.\n";

    Faraday<dimension> faraday{layout, dt};
    Ampere<dimension> ampere{layout};
    Ohm<dimension> ohm{layout};
    Boris<dimension> push{layout, dt};


    auto boundary_condition = BoundaryConditionFactory<dimension>::create("periodic", layout);

    ampere(B, J);
    boundary_condition->fill(J);
    deposit(particles, N, F, *layout);
    boundary_condition->fill(N);
    bulk_velocity<dimension>(N, F, V);
    boundary_condition->fill(V);
    ohm(B, J, N, V, Enew);

    {
        std::string filename = "B.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        auto time_str = to_string_fixed_width(0., 10, 0);
        file.createDataSet("/t/" + time_str + "/B/x", B.x.data());
        file.createDataSet("/t/" + time_str + "/B/y", B.y.data());
        file.createDataSet("/t/" + time_str + "/B/z", B.z.data());
        file.createDataSet("/t/" + time_str + "/J/x", J.x.data());
        file.createDataSet("/t/" + time_str + "/J/y", J.y.data());
        file.createDataSet("/t/" + time_str + "/J/z", J.z.data());
        file.createDataSet("/t/" + time_str + "/E/x", E.x.data());
        file.createDataSet("/t/" + time_str + "/E/y", E.y.data());
        file.createDataSet("/t/" + time_str + "/E/z", E.z.data());
        file.createDataSet("/t/" + time_str + "/V/x", V.x.data());
        file.createDataSet("/t/" + time_str + "/V/y", V.y.data());
        file.createDataSet("/t/" + time_str + "/V/z", V.z.data());
        file.createDataSet("/t/" + time_str + "/N", N.data());
    }

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

        // save particles at t=n before updating to predicted t=n+1
        // because predictor 2 will re-start from t=n
        auto save{particles};
        push(particles, Eavg, Bavg);
        boundary_condition->particles(particles);
        deposit(particles, N, F, *layout);
        boundary_condition->fill(N);
        bulk_velocity<dimension>(N, F, V);
        boundary_condition->fill(V);

        // predictor 2
        faraday(B, Eavg, Bnew);
        ampere(Bnew, J);
        boundary_condition->fill(J);
        ohm(Bnew, J, N, V, Enew);
        std::cout << "post ohm\n";
        average(E, Enew, Eavg);
        average(B, Bnew, Bavg);
        boundary_condition->fill(Eavg);
        boundary_condition->fill(Bavg);

        particles = save;
        push(particles, Eavg, Bavg);
        boundary_condition->particles(particles);
        deposit(particles, N, F, *layout);
        boundary_condition->fill(N);
        bulk_velocity<dimension>(N, F, V);
        boundary_condition->fill(V);

        // corrector
        faraday(B, Eavg, B);
        ampere(B, J);
        boundary_condition->fill(J);
        ohm(B, J, N, V, E);
        boundary_condition->fill(E);


        time += dt;

        if (true)
        {
            std::string filename = "B.h5";
            HighFive::File file(filename, HighFive::File::ReadWrite);
            auto time_str = to_string_fixed_width(time, 10, 0);
            file.createDataSet("/t/" + time_str + "/B/x", B.x.data());
            file.createDataSet("/t/" + time_str + "/B/y", B.y.data());
            file.createDataSet("/t/" + time_str + "/B/z", B.z.data());
            file.createDataSet("/t/" + time_str + "/E/x", E.x.data());
            file.createDataSet("/t/" + time_str + "/E/y", E.y.data());
            file.createDataSet("/t/" + time_str + "/E/z", E.z.data());
            file.createDataSet("/t/" + time_str + "/V/x", V.x.data());
            file.createDataSet("/t/" + time_str + "/V/y", V.y.data());
            file.createDataSet("/t/" + time_str + "/V/z", V.z.data());
            file.createDataSet("/t/" + time_str + "/N", N.data());
        }
    }


    return 0;
}
