#define BOOST_TEST_MODULE "filters_test"

#include <boost/test/included/unit_test.hpp>
#include "object_container.hpp"
#include "filters.hpp"
#include "sphere.hpp"
#include "cylinder.hpp"

template<typename Toc_>
struct collector
{
    void operator()(typename Toc_::iterator i,
            typename Toc_::position_type::value_type dist)
    {
        std::cout << (*i).second << ", " << dist << std::endl;
    }
};


BOOST_AUTO_TEST_CASE(spheres)
{
    typedef double length_type;
    typedef int key_type;
    typedef sphere<length_type> mapped_type;
    typedef object_container<length_type, key_type, mapped_type> oc_type;
    typedef oc_type::position_type pos;
    oc_type oc(1.0, 10);

    oc.insert(std::make_pair(0, mapped_type(pos(0.2, 0.6, 0.4), 0.15)));
    oc.insert(std::make_pair(1, mapped_type(pos(0.2, 0.7, 0.5), 0.05)));
    oc.insert(std::make_pair(2, mapped_type(pos(0.9, 0.1, 0.4), 0.07)));
    oc.insert(std::make_pair(3, mapped_type(pos(0.9, 0.95, 0.4), 0.1)));

    collector<oc_type> col;
    // Construct an iterator for object 1 in oc.
    oc_type::const_iterator f(oc.find(1));
    // Collect all objects from oc that lie within the radius of object 1.
    // The shell of object 0 is at a distance sqrt(2*0.1*0.1) - 0.15 = 
    // -0.0087864..
    // The shell of the objects itself is at a distance -0.05 (Yes it returns 
    // itself also).
    // Objects 2 and 3 should not be found.
    take_neighbor(oc, col, (*f).second);

    // One could check for collisions by asserting that the distance returned 
    // is always > 0.
}

BOOST_AUTO_TEST_CASE(cylinders)
{
    typedef double length_type;
    typedef int key_type;
    typedef cylinder<length_type> mapped_type;
    typedef sphere<length_type> sphere_type;
    typedef object_container<length_type, key_type, mapped_type> oc_type;
    typedef oc_type::position_type pos;
    oc_type oc(1.0, 10);

    oc.insert(std::make_pair(0, mapped_type(pos(0.2, 0.7, 0.7), 0.15, pos(0,0,1), 0.15)));
    oc.insert(std::make_pair(1, mapped_type(pos(0.2, 0.7, 0.5), 0.50, pos(0,0,1), 0.05)));
    oc.insert(std::make_pair(2, mapped_type(pos(0.9, 0.1, 0.4), 0.07, pos(0,0,1), 0.07)));
    oc.insert(std::make_pair(3, mapped_type(pos(0.9, 0.95, 0.4), 0.1, pos(0,0,1), 0.1)));

    cout << std::endl << std::endl;
    typedef collector<oc_type> col_type;
    col_type col;
    // Collect all objects from oc that lie within the radius of object 1.
    //oc_type::const_iterator f(oc.find(1));
    //mapped_type c = (*f).second;
    //take_neighbor<oc_type, col_type, cylinder_type>(oc, col, c);

    // Collect all cylinders from oc who (partly) lie within the radius of 
    // sphere s.
    sphere_type s = sphere_type(pos(0.2, 0.7, 0.5), 0.05);
    take_neighbor<oc_type, col_type, sphere_type>(oc, col, s);
    // Distance to cylinder 0 should be 0.05.
    // Note that the distance to cylinder 1 is not correct because the point 
    // lies inside the cylinder. It should compute the distance to the caps, 
    // but it computes the distance to the radial edge.
}
