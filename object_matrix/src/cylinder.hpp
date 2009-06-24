
#ifndef CYLINDER_HPP
#define CYLINDER_HPP

#include <ostream>
#include "position.hpp"
#include <cmath>
#include <iostream>
//using namespace std;

// Todo. Make sure cylinder is never larger than 1 cellsize or something.  
// Important.
template<typename T_>
struct cylinder
{
    typedef T_ value_type;
    typedef ::position<T_> position_type;
    typedef T_ length_type;

    cylinder()
        : origin(), radius(0), orientation(), size(0) {}

    cylinder(const position_type& _origin, const length_type& _radius, const position_type& _orientation, const length_type& _size )
        : origin(_origin), radius(_radius), orientation(_orientation), size(_size)
    { }

  private:
    length_type calculateDistanceToSelfWithNewOrigin( position_type pos, position_type newOrigin )
    {
	/* First compute the (z,r) components of pos in a coordinate system 
	 * defined by the vectors zVector and rVector, where zVector is the 
	 * orientation of the cylinder (using newOrigin as it's centre) and 
	 * rVector is choosen such that zVector and rVector define a plane in 
	 * which pos lies. */
        position_type posVector = pos - newOrigin;

        length_type z = posVector.dot_product( orientation ); // Can be <0.
        position_type zUnitVector = orientation; // By definition.
        position_type zVector = zUnitVector.scale( z );

        position_type rVector = posVector - zVector;
        length_type r = rVector.length();
        position_type rUnitVector = rVector.scale( 1. / r );


	/* Then compute distance to cylinder. */
        length_type dz = fabs(z) - size;
	length_type distance;
	//cout << "dz = " << dz << endl;
        if(dz > 0){
	    // pos is (either) to the right or to the left of the cylinder.
            if(r < radius){
                distance = dz;
	    }
            else{
		// Difficult case. Compute distance to edge.
                position_type edgeVector = zUnitVector.scale(size) + rUnitVector.scale(radius);
                distance = (posVector - edgeVector).length();
	    }
	}
        else{
	    // pos is somewhere 'parellel' to the cylinder.
	    //cout << " r = " << r;
	    //cout << " radius = " << radius;
            distance = r - radius;
	    //cout << " distance = " << distance;
	}
	//cout << "distance to cylinder = " << distance << endl;
        return distance;
    }

  public:
    length_type calculateDistanceToSelf( position_type pos )
    {
	return calculateDistanceToSelfWithNewOrigin( pos, origin );
    }

    length_type calculateDistanceToSelfWithOffset( position_type pos, position_type offset )
    {
	//cout << "offset " << offset;
	// Because this cylinder is on the other side of one of the periodic 
	// boundaries compared to pos, add offset (for example (-L, 0, 0) to 
	// origin before calculating distance between this cylinder and pos.  
        return calculateDistanceToSelfWithNewOrigin( pos, origin + offset );
    }

    bool operator==(const cylinder& rhs) const
    {
        return origin == rhs.origin && radius == rhs.radius && orientation == rhs.orientation && size == rhs.size;
    }

    bool operator!=(const cylinder& rhs) const
    {
        return !operator==(rhs);
    }

    // Thse need to be public for comparison operator==.
    position_type origin; // centre.
    length_type radius;
    position_type orientation; // should be normalized.
    length_type size; // half length.
};

template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const cylinder<T_>& v)
{
    strm << "{" << v.origin <<  ", " << v.radius << ", " << v.orientation << ", " << v.size << "}";
    return strm;
}

#endif /* CYLINDER_HPP */
