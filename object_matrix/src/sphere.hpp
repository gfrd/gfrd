#ifndef SPHERE_HPP
#define SPHERE_HPP

#include <ostream>
#include "position.hpp"

//// A sphere has a position and a radius.
template<typename T_>
struct sphere
{
    typedef T_ value_type;
    typedef ::position<T_> position_type;
    typedef T_ length_type;

    sphere()
        : position(), radius(0) {}

    //// Initialization list.
    sphere(const position_type& _position, const length_type& _radius)
        : position(_position), radius(_radius) {}

    length_type calculateDistanceToSelf( position_type pos )
    {
        return pos.distance(position) - radius;
    }

    length_type calculateDistanceToSelfWithOffset( position_type pos, position_type offset )
    {
        return pos.distance(position + offset) - radius;
    }

    bool operator==(const sphere& rhs) const
    {
        return position == rhs.position && radius == rhs.radius;
    }

    bool operator!=(const sphere& rhs) const
    {
        return !operator==(rhs);
    }
    
    position_type position;
    length_type radius;
};

template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const sphere<T_>& v)
{
    strm << "{" << v.position <<  ", " << v.radius << "}";
    return strm;
}

#endif /* SPHERE_HPP */
