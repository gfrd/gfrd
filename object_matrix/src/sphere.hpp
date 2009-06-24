#ifndef SPHERE_HPP
#define SPHERE_HPP

#include <ostream>
#include "position.hpp"

/* A sphere has an origin and a size (radius). */
template<typename T_>
struct sphere
{
    typedef T_ value_type;
    typedef ::position<T_> position_type;
    typedef T_ length_type;

    sphere()
        : origin(), size(0) {}

    sphere(const position_type& _origin, const length_type& _size)
        : origin(_origin), size(_size) {}

    length_type calculateDistanceToSelf( position_type pos )
    {
        return pos.distance(origin) - size;
    }

    length_type calculateDistanceToSelfWithOffset( position_type pos, position_type offset )
    {
	// Because this sphere is on the other side of one of the periodic 
	// boundaries compared to pos, add offset (for example (-L, 0, 0) to 
	// origin before calculating distance between this sphere and pos. 
        return pos.distance(origin + offset) - size;
    }

    bool operator==(const sphere& rhs) const
    {
        return origin == rhs.origin && size == rhs.size;
    }

    bool operator!=(const sphere& rhs) const
    {
        return !operator==(rhs);
    }
    
    position_type origin;
    length_type size;
};

template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const sphere<T_>& v)
{
    strm << "{" << v.origin <<  ", " << v.size << "}";
    return strm;
}

#endif /* SPHERE_HPP */
