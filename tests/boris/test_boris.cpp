#include "pusher.hpp"

#include "highfive/highfive.hpp"

#include <iostream>
#include <vector>

void uniform_bz()
{
    Particle<1> particle;
    particle.position[0] = 5.05;
    particle.v[0]        = 0.0;
    particle.v[1]        = 2.0;
    particle.v[2]        = 0.0;
    particle.weight      = 1.0;
    particle.mass        = 1.0;
    particle.charge      = 1.0;
    std::vector<Particle<1>> particles{particle};

    double time                     = 0.;
    double final_time               = 3.141592 * 4;
    double dt                       = 0.001;
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {1000};
    std::array<double, dimension> cell_size      = {0.1};
    auto constexpr nbr_ghosts                    = 1;
    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> E{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};


    for (auto ix = layout->dual_dom_start(Direction::X); ix <= layout->dual_dom_end(Direction::X);
         ++ix)
    {
        B.z(ix) = 3.0; // Uniform magnetic field in z-direction
    }


    // Create a Boris pusher
    Boris<dimension> push{layout, dt};
    std::vector<double> x;
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;


    while (time < final_time)
    {
        x.push_back(particles[0].position[0]);
        vx.push_back(particles[0].v[0]);
        vy.push_back(particles[0].v[1]);
        vz.push_back(particles[0].v[2]);
        // Push the particles using the Boris pusher
        push(particles, E, B);


        double const iCell_float = particle.position[0] / layout->cell_size(Direction::X)
                                   + layout->dual_dom_start(Direction::X);
        int const iCell = static_cast<int>(iCell_float);
        if (iCell > layout->dual_dom_end(Direction::X))
        {
            particle.position[0] -= layout->dom_size(Direction::X);
        }
        else if (iCell < layout->dual_dom_end(Direction::X))
        {
            particle.position[0] += layout->dom_size(Direction::X);
        }

        time += dt;
    }

    {
        std::string filename = "uniform_bz.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/x", x);
        file.createDataSet("/vx", vx);
        file.createDataSet("/vy", vy);
        file.createDataSet("/vz", vz);
    }

    std::cout << "This is a test for the Boris pusher in a hybrid simulation.\n";
}

void drift_ey()
{
    Particle<1> particle;
    particle.position[0] = 5.05;
    particle.v[0]        = 0.0;
    particle.v[1]        = 1.0;
    particle.v[2]        = 0.0;
    particle.weight      = 1.0;
    particle.mass        = 1.0;
    particle.charge      = 1.0;
    std::vector<Particle<1>> particles{particle};

    double time                     = 0.;
    double final_time               = 3.141592 * 4;
    double dt                       = 0.001;
    std::size_t constexpr dimension = 1;

    std::array<std::size_t, dimension> grid_size = {1000};
    std::array<double, dimension> cell_size      = {0.1};
    auto constexpr nbr_ghosts                    = 1;
    auto layout = std::make_shared<GridLayout<dimension>>(grid_size, cell_size, nbr_ghosts);

    VecField<dimension> E{layout, {Quantity::Ex, Quantity::Ey, Quantity::Ez}};
    VecField<dimension> B{layout, {Quantity::Bx, Quantity::By, Quantity::Bz}};


    for (auto ix = layout->dual_dom_start(Direction::X); ix <= layout->dual_dom_end(Direction::X);
         ++ix)
    {
        B.z(ix) = 8.;  // Uniform magnetic field in z-direction
        E.y(ix) = 1.0; // Uniform magnetic field in z-direction
    }


    // Create a Boris pusher
    Boris<dimension> push{layout, dt};
    std::vector<double> x;
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;


    while (time < final_time)
    {
        // Push the particles using the Boris pusher
        push(particles, E, B);

        x.push_back(particles[0].position[0]);
        vx.push_back(particles[0].v[0]);
        vy.push_back(particles[0].v[1]);
        vz.push_back(particles[0].v[2]);

        double const iCell_float = particle.position[0] / layout->cell_size(Direction::X)
                                   + layout->dual_dom_start(Direction::X);
        int const iCell = static_cast<int>(iCell_float);
        if (iCell > layout->dual_dom_end(Direction::X))
        {
            particle.position[0] -= layout->dom_size(Direction::X);
        }
        else if (iCell < layout->dual_dom_end(Direction::X))
        {
            particle.position[0] += layout->dom_size(Direction::X);
        }

        time += dt;
    }

    {
        std::string filename = "drift_ey.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/x", x);
        file.createDataSet("/vx", vx);
        file.createDataSet("/vy", vy);
        file.createDataSet("/vz", vz);
    }

    std::cout << "This is a test for the Boris pusher in a hybrid simulation.\n";
}

int main()
{
    uniform_bz();
    drift_ey();
}
