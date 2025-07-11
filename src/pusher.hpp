#ifndef PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP
#define PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP


#include "particle.hpp"


#include <cstddef>
#include <vector>

namespace HYBIRT
{
template<std::size_t dim>
class Pusher
{
protected:
public:
    virtual void move(std::vector<Particle<dim>>& particles) = 0;

    virtual ~Pusher() {}
};

} // namespace HYBIRT

#endif
