
#ifndef CYLINDER_HPP
#define CYLINDER_HPP

#include <ostream>
#include "position.hpp"
#include <cmath>

// Todo. Make sure cylinder is never larger than 1 cellsize or something.  
template<typename T_>
struct cylinder
{
    typedef T_ value_type;
    typedef ::position<T_> position_type;
    typedef T_ length_type;

    cylinder()
        : origin(), radius(0), orientationZ(), size(0) {}

    cylinder(const position_type& _origin, const length_type& _radius, const position_type& _orientationZ, const length_type& _size )
        : origin(_origin), radius(_radius), orientationZ(_orientationZ), size(_size)
    { }

  private:
    length_type calculateDistanceToSelfWithNewOrigin( position_type pos, position_type newOrigin )
    {
	/* First compute the (z,r) components of pos in a coordinate system 
	 * defined by the vectors unitR and orientationZ, where unitR is
	 * choosen such that unitR and orientationZ define a plane in which
	 * pos lies. */
        position_type posVector = pos - newOrigin;

        length_type z = posVector.dot_product( orientationZ ); // Can be <0.
        position_type posVectorZ = orientationZ.scale( z );

        position_type posVectorR = posVector - posVectorZ;
        length_type r = posVectorR.length();


	/* Then compute distance to cylinder. */
        length_type dz = fabs(z) - size;
	length_type dr = r - radius;
	length_type distance;
        if(dz > 0){
	    // pos is (either) to the right or to the left of the cylinder.
            if(r > radius){
		// Compute distance to edge.
		distance = std::sqrt( dz*dz + dr*dr );
	    }
            else{
                distance = dz;
	    }
	}
        else{
	    if(dr > radius){
		// pos is somewhere 'parellel' to the cylinder.
		distance = dr;
	    }
	    else{
		// Inside cylinder. 
		distance = std::max(dr, dz);
	    }
	}
        return distance;
    }

  public:
    length_type calculateDistanceToSelf( position_type pos )
    {
	return calculateDistanceToSelfWithNewOrigin( pos, origin );
    }

    length_type calculateDistanceToSelfWithOffset( position_type pos, position_type offset )
    {
	// Because this cylinder is on the other side of one of the periodic 
	// boundaries compared to pos, add offset (for example (-L, 0, 0) to 
	// origin before calculating distance between this cylinder and pos.  
        return calculateDistanceToSelfWithNewOrigin( pos, origin + offset );
    }

    bool operator==(const cylinder& rhs) const
    {
        return origin == rhs.origin && radius == rhs.radius && orientationZ == rhs.orientationZ && size == rhs.size;
    }

    bool operator!=(const cylinder& rhs) const
    {
        return !operator==(rhs);
    }

    // Thse need to be public for comparison operator==.
    position_type origin; // centre.
    length_type radius;
    position_type orientationZ; // should be normalized.
    length_type size; // half length.
};

template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const cylinder<T_>& v)
{
    strm << "{" << v.origin <<  ", " << v.radius << ", " << v.orientationZ << ", " << v.size << "}";
    return strm;
}

#endif /* CYLINDER_HPP */
