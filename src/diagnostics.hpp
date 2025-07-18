#ifndef HYBIRT_DIAGNOSTICS_HPP
#define HYBIRT_DIAGNOSTICS_HPP

#include "field.hpp"
#include "vecfield.hpp"
#include "particle.hpp"
#include "population.hpp"

#include "highfive/highfive.hpp"

#include <iomanip>
#include <vector>
#include <string>

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




template<std::size_t dim>
void diags_write_fields(VecField<dim> const& B, VecField<dim> const& E, VecField<dim> const& V,
                        Field<dim> const& N, double time,
                        HighFive::File::AccessMode mode = HighFive::File::ReadWrite)
{
    std::string filename = "fields.h5";
    HighFive::File file(filename, mode);
    auto const time_str = to_string_fixed_width(time, 10, 0);
    file.createDataSet("/t/" + time_str + "/Bx", B.x.data());
    file.createDataSet("/t/" + time_str + "/By", B.y.data());
    file.createDataSet("/t/" + time_str + "/Bz", B.z.data());
    file.createDataSet("/t/" + time_str + "/Ex", E.x.data());
    file.createDataSet("/t/" + time_str + "/Ey", E.y.data());
    file.createDataSet("/t/" + time_str + "/Ez", E.z.data());
    file.createDataSet("/t/" + time_str + "/Vx", V.x.data());
    file.createDataSet("/t/" + time_str + "/Vy", V.y.data());
    file.createDataSet("/t/" + time_str + "/Vz", V.z.data());
    file.createDataSet("/t/" + time_str + "/N", N.data());
}



template<std::size_t dim>
void diags_write_particles(std::vector<Population<dim>> const& populations, double time,
                           HighFive::File::AccessMode mode = HighFive::File::ReadWrite)
{
    for (auto const& pop : populations)
    {
        std::string filename = "particles_" + pop.name() + ".h5";
        HighFive::File file(filename, mode);

        std::vector<double> x, y, z, vx, vy, vz;
        x.reserve(pop.particles().size());
        vx.reserve(pop.particles().size());
        vy.reserve(pop.particles().size());
        vz.reserve(pop.particles().size());
        for (auto const& particle : pop.particles())
        {
            x.push_back(particle.position[0]);
            vx.push_back(particle.v[0]);
            vy.push_back(particle.v[1]);
            vz.push_back(particle.v[2]);
        }
        auto const time_str = to_string_fixed_width(time, 10, 0);
        file.createDataSet("/t/" + time_str + "/x", x);
        file.createDataSet("/t/" + time_str + "/vx", vx);
        file.createDataSet("/t/" + time_str + "/vy", vy);
        file.createDataSet("/t/" + time_str + "/vz", vz);
    }
}


#endif
