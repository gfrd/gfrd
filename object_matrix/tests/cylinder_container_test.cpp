#define BOOST_TEST_MODULE "object_container_test"

#include <functional>
#include <iostream>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "object_container.hpp"
#include "cylinder.hpp"

BOOST_AUTO_TEST_CASE(insert)
{
    typedef double length_type;
    typedef int key_type;
    typedef cylinder<length_type> mapped_type;
    typedef object_container<length_type, key_type, mapped_type> oc_type;
    typedef oc_type::position_type pos;
    oc_type oc(1000.0, 10);
    BOOST_CHECK_CLOSE(100., oc.cell_size(), 0.001);
    {
        std::pair<oc_type::iterator, bool> ir(
                oc.insert(std::make_pair(
                    0, oc_type::mapped_type(pos(500, 500, 500), 25, pos(0, 0, 1), 1000))));
        BOOST_CHECK_EQUAL(true, ir.second);  // Normal insert.
        BOOST_CHECK(oc.end() != oc.find(0)); // Key 0 exists.
        BOOST_CHECK(oc.end() == oc.find(1)); // Key 1 doesn't exist.
    }
    {
	// Update.
        std::pair<oc_type::iterator, bool> ir(
                oc.insert(std::make_pair(
                    0, oc_type::mapped_type(pos(500, 500, 500), 25, pos(0, 0, 1), 1000))));
        BOOST_CHECK_EQUAL(false, ir.second); // False: this was an update.
	// ir.first is an iterator to the value you inserted. So accessing 
	// it's second element should return you the object.
        BOOST_CHECK_EQUAL(oc_type::mapped_type(pos(500, 500, 500), 25, pos(0, 0, 1), 1000),
                (*ir.first).second);
        BOOST_CHECK(oc.end() != oc.find(0));
        BOOST_CHECK(oc.end() == oc.find(1));
    }
    {
	// Another update.
        std::pair<oc_type::iterator, bool> ir(
                oc.insert(std::make_pair(
                    0, oc_type::mapped_type(pos(500, 500, 500), 25, pos(0, 0, 1), 1000))));
        BOOST_CHECK_EQUAL(false, ir.second);
        BOOST_CHECK_EQUAL(oc_type::mapped_type(pos(500, 500, 500), 25, pos(0, 0, 1), 1000),
                (*ir.first).second);
        BOOST_CHECK(oc.end() != oc.find(0));
        BOOST_CHECK(oc.end() == oc.find(1));
    }
}

template<typename Toc_>
struct collector
{
    void operator()(typename Toc_::iterator i)
    {
        std::cout << (*i).second << std::endl;
    }
};

template<typename Toc_>
struct collector2
{
    void operator()(typename Toc_::iterator i,
            const typename Toc_::position_type& pos_off)
    {
        std::cout << (*i).second;
	// Signal if periodic boundary condition is applied (i.e. (-1,0,0) 
	// etc). 
        std::cout << pos_off << std::endl;
    }
};

BOOST_AUTO_TEST_CASE(each_neighbor)
{
    typedef double length_type;
    typedef int key_type;
    typedef cylinder<length_type> mapped_type;
    typedef object_container<length_type, key_type, mapped_type> oc_type;
    typedef oc_type::position_type pos;
    oc_type oc(1000., 10);
    BOOST_CHECK_CLOSE(100., oc.cell_size(), 0.001);

    // Insert 4 values.
    // Range of cells: [0, 0.1), [0.1, 0.2), .., [0.8, 0.9)
    // [0.9, 1) = [0.9, 0)
    oc.insert(std::make_pair(0, oc_type::mapped_type(pos(500, 500, 0), 25, pos(0,0,1), 50)));
    BOOST_CHECK(oc.end() != oc.find(0));
    BOOST_CHECK(oc.end() == oc.find(1));

    collector<oc_type> col;
    // Should return values 1.
    oc.each_neighbor(oc.index(pos(500, 500, 100)), col);
    std::cout << "--" << std::endl;

    // No periodic boundary condition. Should return no values.
    oc.each_neighbor(oc.index(pos(500, 500, 1000)), col);
    std::cout << "--" << std::endl;

    collector2<oc_type> col2;
    // Periodic boundary condition. Should return element 1 after applying 
    // periodic boundary condition in z (add 1000 to z coordinate of the 
    // origin of the cylinder to be in the same neighbourhood as reference 
    // point), so: (0,0,1000).
    oc.each_neighbor_cyclic(oc.index(pos(500, 500, 900)), col2);
    std::cout << "--" << std::endl;

    oc.insert(std::make_pair(1, oc_type::mapped_type(pos(500, 500, 900), 25, pos(0,0,1), 50)));
    BOOST_CHECK(oc.end() != oc.find(0));
    BOOST_CHECK(oc.end() != oc.find(1));
    BOOST_CHECK(oc.end() == oc.find(2));

    // Periodic boundary condition. Should return element 0 (0, 0, 0) and 
    // element 1 (0,0,-1000).
    oc.each_neighbor_cyclic(oc.index(pos(500, 500, 0)), col2);
}

