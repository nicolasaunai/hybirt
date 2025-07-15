#ifndef HYBRIDIR_FIELD_HPP
#define HYBRIDIR_FIELD_HPP

#include <cstddef>
#include <vector>
#include <numeric>

template<std::size_t dimension>
class Field
{
public:
    Field(std::array<std::size_t, dimension> grid_size)
        : m_size{grid_size}
        , m_data(std::accumulate(grid_size.begin(), grid_size.end(), 1,
                                 std::multiplies<std::size_t>()),
                 1.0)
    {
    }


    template<typename... Indexes>
    double& operator()(Indexes... ijk)
    {
        auto indexes = std::forward_as_tuple(ijk...);

        if constexpr (dimension == 1)
        {
            auto index = std::get<0>(indexes);
            return m_data[index];
        }
    }

    template<typename... Indexes>
    double const& operator()(Indexes... ijk) const
    {
        auto indexes = std::forward_as_tuple(ijk...);

        if constexpr (dimension == 1)
        {
            auto index = std::get<0>(indexes);
            return m_data[index];
        }
    }

private:
    std::array<std::size_t, dimension> m_size;
    std::vector<double> m_data;
};


#endif // HYBRIDIR_VECFIELD_HPP
