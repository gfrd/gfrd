#ifndef POSITION_HPP
#define POSITION_HPP

#include <ostream>
#include <functional>
#include <algorithm>
#include <cmath>

//// An array class that has the same interface as STL containers, but unlike
//// vector has a fixed number of elements (so faster).
#include <boost/array.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>


//// A position is an array of lenght 3 of type T_. 
template<typename T_>
struct position: public boost::array<T_, 3>
{
    typedef boost::array<T_, 3> base_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::size_type size_type;

    position()
    {
        //// This works because position derives from boost::array.
        (*this)[0] = 0;
        (*this)[1] = 0;
        (*this)[2] = 0;
    }

    //// constructors
    //// 1. If a is a reference to an array of type T_, cast it to a pointer 
    //// to a base_type (why the extra *?), which is a boost array of length 3  
    //// in this case, and instantiate such an array (initialization list).
    position(const T_ (&a)[3]): base_type(
            *reinterpret_cast<const base_type*>(&a)) {}

    //// 2. argument is an array of type T_
    position(const T_ a[3]): base_type(
            *reinterpret_cast<const base_type*>(a)) {}

    //// 3. argument is a reference to a boost array
    position(const base_type& a): base_type(a) {}

    //// 4. 3 value arguments specified.
    position(value_type x, value_type y, value_type z)
    {
        (*this)[0] = x;
        (*this)[1] = y;
        (*this)[2] = z;
    }

    //// Return a reference to the x value (instead of returning a copy, 
    //// faster).
    value_type& x()
    {
        return (*this)[0];
    }

    const value_type& x() const
    {
        return (*this)[0];
    }

    value_type& y()
    {
        return (*this)[1];
    }

    const value_type& y() const
    {
        return (*this)[1];
    }

    value_type& z()
    {
        return (*this)[2];
    }

    const value_type& z() const
    {
        return (*this)[2];
    }


    const value_type distance_sq(const position& that) const
    {
        return std::pow((*this)[0] - that[0], 2)
            + std::pow((*this)[1] - that[1], 2)
            + std::pow((*this)[2] - that[2], 2);
    }

    const value_type distance(const position& that) const
    {
        return std::sqrt(distance_sq(that));
    }

    // Todo. Why are we defining our own vector class really? Should use some 
    // library.
    value_type dot_product(const position& that) const
    {
        position temp;
        // Take element-wise product.
        std::transform(
            boost::const_begin(*this), boost::const_end(*this),
            boost::const_begin(that), boost::begin(temp),
            std::multiplies<value_type>());

        // Sum the result.
        value_type sum = 0;
        for(typename base_type::iterator it = temp.begin(); it != temp.end(); ++it) {
           sum += (*it); 
        }

        return sum;
    }

    // Can someone explain to me why we don't just do:
    // retval[0] = (*this)[0] + (*that)[0];
    // retval[1] = (*this)[1] + (*that)[1];
    // retval[2] = (*this)[2] + (*that)[2];
    position operator+(const position& that) const
    {
        position retval;
        std::transform(
            boost::const_begin(*this), boost::const_end(*this),
            boost::const_begin(that), boost::begin(retval),
            std::plus<value_type>());
        return retval;
    }

    position operator-(const position& that) const
    {
        position retval;
        std::transform(
            boost::const_begin(*this), boost::const_end(*this),
            boost::const_begin(that), boost::begin(retval),
            std::minus<value_type>());
        return retval;
    }

    position scale(const value_type multiplier) const
    {
        position retval;
        retval[0] = (*this)[0] * multiplier;
        retval[1] = (*this)[1] * multiplier;
        retval[2] = (*this)[2] * multiplier;
        return retval;
    }
};

template<typename Tstrm_, typename T_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const position<T_>& v)
{
    strm << "(" << v.x() <<  ", " << v.y() <<  ", " << v.z() << ")";
    return strm;
}

#endif /* POSITION_HPP */
