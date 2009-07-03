
#ifndef BOX_HPP
#define BOX_HPP

#include <ostream>
#include "position.hpp"
#include <cmath>

// Todo. Make sure box is never larger than 1 cellsize or something.  
template<typename T_>
struct box
{
    typedef T_ value_type;
    typedef ::position<T_> position_type;
    typedef T_ length_type;

    box()
        : origin(), unitX(), unitY(), unitZ(), Lx(), Ly(), Lz() {}

    box(const position_type& _origin, const position_type& _unitX, const position_type& _unitY, const position_type& _unitZ, const length_type _Lx, const length_type _Ly, const length_type _Lz )
        : origin(_origin), unitX(_unitX), unitY(_unitY), unitZ(_unitZ), Lx(_Lx), Ly(_Ly), Lz(_Lz)
    {
    }

  private:
    length_type calculateDistanceToSelfWithNewOrigin( position_type pos, position_type newOrigin )
    {
	/* First compute the (x,y,z) components of pos in a coordinate system 
	 * defined by the vectors unitX, unitY, unitZ,
	 * where. */
        position_type posVector = pos - newOrigin;

        length_type x = posVector.dot_product( unitX ); // Can be <0.
        length_type y = posVector.dot_product( unitY ); // Can be <0.
        length_type z = posVector.dot_product( unitZ ); // Can be <0.

	/* Then compute distance to box. */
        length_type dx = fabs(x) - Lx;
        length_type dy = fabs(y) - Ly;
        length_type dz = fabs(z) - Lz;

	length_type distance;

        if(dx > 0){
            if(dy > 0){
		if(dz > 0){
		    distance = std::sqrt( dx*dx + dy*dy + dz*dz );
		}
		else{
		    distance = std::sqrt( dx*dx + dy*dy );
		}
	    }
	    else{
		if(dz > 0){
		    distance = std::sqrt( dx*dx + dz*dz );
		}
		else{
		    distance = dx;
		}
	    }
	}
	else{
            if(dy > 0){
		if(dz > 0){
		    distance = std::sqrt( dy*dy + dz*dz );
		}
		else{
		    distance = dy;
		}
	    }
	    else{
		if(dz > 0){
		    distance = dz;
		}
		else{
		    // Inside box. Pick negative distance closest to 0.
		    distance = std::max(std::max(dx,dy), dz);
		}
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
	// Because this box is on the other side of one of the periodic 
	// boundaries compared to pos, add offset (for example (-L, 0, 0) to 
	// origin before calculating distance between this box and pos.  
        return calculateDistanceToSelfWithNewOrigin( pos, origin + offset );
    }

    bool operator==(const box& rhs) const
    {
	return origin == rhs.origin && unitX == rhs.unitX &&
	    unitY == rhs.unitY && unitZ == rhs.unitZ &&
	    Lx == rhs.Lx && Ly == rhs.Ly  && Lz == rhs.Lz;
    }

    bool operator!=(const box& rhs) const
    {
        return !operator==(rhs);
    }

    position_type origin; // centre.
    // Unit orientation vectors.
    position_type unitX;
    position_type unitY;
    position_type unitZ;
    length_type Lx;
    length_type Ly;
    length_type Lz;
};

template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const box<T_>& v)
{
    strm << "{" << v.origin <<  ", " << v.Lx*v.unitX << ", " << v.Ly*v.unitY << ", " << v.Lz*v.unitZ << "}";
    return strm;
}

#endif /* BOX_HPP */
