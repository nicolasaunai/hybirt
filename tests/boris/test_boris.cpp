#include "pusher.hpp"

#include "highfive/highfive.hpp"

#include <iostream>
#include <vector>

int main()
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
        B.z(ix) = 1.0; // Uniform magnetic field in z-direction
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

        time += dt;
    }

    {
        std::string filename = "orbit.h5";
        HighFive::File file(filename, HighFive::File::Truncate);
        file.createDataSet("/x", x);
        file.createDataSet("/vx", vx);
        file.createDataSet("/vy", vy);
        file.createDataSet("/vz", vz);
    }

    std::cout << "This is a test for the Boris pusher in a hybrid simulation.\n";
    return 0;
}
