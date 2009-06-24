#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include <functional>
#include <cmath>
#include "position.hpp"
#include "sphere.hpp"
#include <iostream>
//#include <typeinfo>
/* See cyclic_neighbor_filter below for explanation.
 * This does the same but without applying periodic boundary conditions. */
template<typename Toc_, typename Tfun_, typename Tsphere_>
class neighbor_filter
        : public std::unary_function<typename Toc_::reference, void>
{
    typedef typename Toc_::iterator argument_type;
    typedef void result_type;

public:
    inline neighbor_filter(Tfun_& next, const Tsphere_& cmp)
        : next_(next), cmp_(cmp) {}

    inline result_type operator()(argument_type i) const
    {
        typename argument_type::reference item(*i);

	/* Todo. This doesn't work anymore if mapped_type is a cylinder.
        if (cmp_ == item.second)
        {
            return;
        }
	*/

        typename Toc_::mapped_type object = item.second;
        const double dist(
            // FIXME: something's wrong
	    //const_cast<position<double>& > 
            (object).calculateDistanceToSelf(cmp_.origin));
        if (dist < cmp_.size)
        {
            next_(i, dist);
        }
    }

private:
    Tfun_& next_;
    const Tsphere_& cmp_;
};

template<typename Toc_, typename Tfun_, typename Tsphere_>
inline void take_neighbor(Toc_& oc, Tfun_& fun,
        Tsphere_& cmp)
{
    //std::cout << std::endl << typeid( oc ).name() << std::endl << std::endl; 
    //std::cout << typeid( fun ).name() << std::endl << std::endl;
    //std::cout << typeid( cmp ).name() << std::endl << std::endl; 
    oc.each_neighbor(oc.index(cmp.origin),
                     neighbor_filter<Toc_, Tfun_, Tsphere_>(fun, cmp));
}

// Todo. Here we assume cmp is a sphere!
template<typename Toc_, typename Tfun_, typename Tsphere_>
class cyclic_neighbor_filter
        : public std::binary_function<
            typename Toc_::reference,
            const typename Toc_::position_type&,
            void>
{
    typedef typename Toc_::iterator first_argument_type;
    typedef const typename Toc_::position_type& second_argument_type;
    typedef void result_type;

public:
    // Called via some steps from neighbors_array().
    // Arguments are:
    //	+ next=original collector in which results will be stored.
    //	+ cmp=sphere with position and radius.
    // Objects within the sphere are selected.
    inline cyclic_neighbor_filter(Tfun_& next,
            const Tsphere_& cmp)
        : next_(next), cmp_(cmp) {}

    // Called from each_neighbor_cyclic_loops in object_container.hpp.
    // Arguments are:
    //	+ i=iterator
    //	+ p=position_offset.
    inline result_type operator()(first_argument_type i,
            second_argument_type p) const
    {
        typename first_argument_type::reference item(*i);
        typename Toc_::mapped_type object = item.second;
	
	/* Todo. This doesn't work anymore if mapped_type is a cylinder.
        if (cmp_ == object)
        {
	    // In case cmp_, the 'virtual' sphere', is actually an object in  
	    // our matrix (i.e the matrix also contains a sphere at the same 
	    // position and with the same radius as cmp), don't include it.  
            return;
        }
	*/

	// Calculate distance from position of sphere to the shell of the ith 
	// item.
        const double dist(
            // FIXME: something's wrong
	    //const_cast<position<double>& >(cmp_.position)
            (object).calculateDistanceToSelfWithOffset(cmp_.origin, p));
        if (dist < cmp_.size)
        {
	    // If distance is within the radius of the sphere, add i to 
	    // collector.
            next_(i, dist);
        }
	    // Else: object i is not within the radius of sphere cmp_.
    }

private:
    Tfun_& next_;
    const Tsphere_& cmp_;
};

// Called via some steps from  
// subSpaceSimulator.getNeighborsWithinRadiusNoSort().
template<typename Toc_, typename Tfun_, typename Tsphere_>
inline void take_neighbor_cyclic(Toc_& oc, Tfun_& fun,
         const Tsphere_& cmp)
{
    oc.each_neighbor_cyclic(oc.index(cmp.origin),
            cyclic_neighbor_filter<Toc_, Tfun_, Tsphere_>(fun, cmp));
}

#endif /* ALGORITHM_HPP */
