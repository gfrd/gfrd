
#ifndef CYLINDER_HPP
#define CYLINDER_HPP

#include <ostream>
#include "position.hpp"
#include <cmath>

// Todo.
// Make sure cylinder is never larger than 1 cellsize or something. Important.
template<typename T_>
struct cylinder
{
    typedef T_ value_type;
    typedef ::position<T_> position_type;
    typedef T_ length_type;

    cylinder()
        : position(), orientation(), radius(0), halfLength(0) {}

    cylinder(const position_type& _position, const position_type& _orientation, const length_type& _radius, const length_type& _halfLength )
        : position(_position), orientation(_orientation), radius(_radius), halfLength(_halfLength)
    {
        orientationVector = orientation.scale(halfLength); 
        leftEnd = position - orientationVector;
        rightEnd = position + orientationVector;

    }

    length_type calculateDistanceToSelf( position_type pos )
    {
        // Calculate position relative to origin of cylinder;
        position_type posRelative = pos - position;
        // Calculate inner product of posRelative with orientation of 
        // cylinder.
        length_type dz = posRelative.dot_product( orientation );

        // Todo. Don't use caps, but calculate for real.
        // Assume origin of cylinder is in the centre -> symmetry.
        if( dz < -halfLength ){
            // pos is to the left of cylinder. Place cap and compute distance.
            return posRelative.distance(leftEnd) - radius;

        }
        else if( dz > halfLength ){
            // pos is to the right of cylinder. Place cap and compute 
            // distance.
            return posRelative.distance(rightEnd) - radius;
        }
        else{
            // pos is not to left/right of cylinder.
            position_type posRelativeZvector = orientation.scale(dz); 
            position_type posRelativeRvector = posRelative - posRelativeZvector;
            return posRelativeRvector.distance(position) - radius;
        }
    }

    length_type calculateDistanceToSelfWithOffset( position_type pos, position_type offset )
    {
        // Todo.
        return pos.distance(position + offset) - radius;
    }

    bool operator==(const cylinder& rhs) const
    {
        return position == rhs.position && orientation == rhs.orientation && radius == rhs.radius && halfLength == rhs.halfLength;
    }

    bool operator!=(const cylinder& rhs) const
    {
        return !operator==(rhs);
    }
    
    position_type position;
    position_type orientation;
    position_type orientationVector;
    position_type leftEnd;
    position_type rightEnd;
    length_type radius;
    length_type halfLength;
};

template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const cylinder<T_>& v)
{
    strm << "{" << v.position <<  ", " << v.orientation << ", " << v.radius << ", " << v.halfLength << "}";
    return strm;
}

#endif /* CYLINDER_HPP */
