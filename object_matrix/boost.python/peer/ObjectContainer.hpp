#ifndef OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP
#define OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP



#include "../../../config.h"

#include <functional>
#include <string>
#include <vector>

#if HAVE_UNORDERED_MAP
#include <unordered_map>
#elif HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
#elif HAVE_EXT_HASH_MAP
#include <ext/hash_map>
#else
#include <map>
#endif /* HAVE_UNORDERED_MAP */

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits.hpp>
#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <numpy/arrayobject.h>

#include "peer/utils.hpp"
#include "peer/numpy/type_mappings.hpp"

#include "peer/tuple_converters.hpp"
#include "peer/numpy/pyarray_backed_allocator.hpp"
#include "peer/numpy/ndarray_converters.hpp"

#include "position.hpp"
#include "filters.hpp"
#include "sphere.hpp"
#include "cylinder.hpp"
#include "object_container.hpp"

#if OBJECTMATRIX_USE_ITERATOR
#include <boost/coroutine/generator.hpp>
#include <boost/bind.hpp>
#include "peer/generator_support.hpp"
#endif /* OBJECTMATRIX_USE_ITERATOR */

namespace peer {

//// The definition of the mapper.
template<typename Tkey_, typename Tval_>
struct get_mapper_mf
{
#if HAVE_UNORDERED_MAP
    typedef std::unordered_map<Tkey_, Tval_> type;
#elif HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<Tkey_, Tval_> type;
#elif HAVE_EXT_HASH_MAP
    typedef __gnu_cxx::hash_map<Tkey_, Tval_> type;
#else 
    typedef std::map<Tkey_, Tval_> type;
#endif
};

//// Other definition of the mapper, with (partial) template specialization.  
//// This one is used if key is a python object instead of of type Tkey_.
template<typename Tval_>
struct get_mapper_mf<boost::python::object, Tval_>
{
#if HAVE_UNORDERED_MAP || HAVE_TR1_UNORDERED_MAP || HAVE_EXT_HASH_MAP
    struct hasher: public std::unary_function<boost::python::object, std::size_t>
    {
        typedef boost::python::object argument_type;
        typedef std::size_t result_type;

        result_type operator()(const argument_type& arg) const
        {
            return static_cast<result_type>((long)PyObject_Hash(arg.ptr()));
        }
    };
#endif


#if HAVE_UNORDERED_MAP
    typedef std::unordered_map<boost::python::object, Tval_, hasher> type;
#elif HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<boost::python::object, Tval_, hasher> type;
#elif HAVE_EXT_HASH_MAP
    typedef __gnu_cxx::hash_map<boost::python::object, Tval_, hasher> type;
#else 
    typedef std::map<boost::python::object, Tval_> type;
#endif
};

template<typename Tmapped_type_>
class ObjectContainer
{
public:
    typedef boost::python::object key_type;
    // Mapped_type now a ObjectContainer template argument.
    typedef Tmapped_type_ mapped_type;
    typedef typename mapped_type::length_type length_type;
    typedef typename mapped_type::position_type position_type;
    //// An impl_type instance is an implemented object_container.
    ////    + length_type = usually double
    ////    + key_type = python::object
    ////    + mapped_type = sphere or cylinder
    ////    + mapper = get_mapper_mf
    typedef ::object_container<length_type, key_type, mapped_type, get_mapper_mf> impl_type;
    typedef typename impl_type::size_type size_type;
    typedef typename impl_type::matrix_type::size_type matrix_size_type;

    // sphere_type is only used for finding objects within a certain radius of 
    // a point. It has nothing to do with the used mapped_type.
    typedef const sphere<length_type> sphere_type;

#ifdef OBJECTMATRIX_USE_ITERATOR
    class Generators
    {
    public:
        typedef std::pair<impl_type::iterator, position_type::value_type>
                result_type;

        typedef boost::coroutines::generator<const result_type*> generator_type; 
    private:
        struct collector: public std::binary_function<
                impl_type::reference, position_type::value_type, void>
        {
        public:
            //// These 2 typedefs are not used.
            typedef impl_type::iterator first_argument_type;
            typedef position_type::value_type second_argument_type;
            typedef void result_type;

        public:
            inline collector(generator_type::self& self,
                    impl_type::iterator end)
                    : self_(self), last_(end, 0)
            {
            }

            inline void operator()(impl_type::iterator i,
                    const position_type::value_type& d)
            {
                last_.first = i;
                last_.second = d;
                self_.yield(&last_);
            }

        private:
            generator_type::self& self_;
            Generators::result_type last_;
        };

    public:
        //// Todo.
        //// When converting from spheres to cylinders, don't care about these  
        //// 2 methods, they are called from a part we're not using.
        inline static Enumerator<const result_type&>* enumerate_neighbors(
            ObjectContainer<sphere<double>>& cntnr, const typename mapped_type& sphere)
        {
            return util::make_enumerator(generator_type(
                boost::bind(&gen_take_neighbor, _1, cntnr, sphere)));
        }

        inline static Enumerator<const result_type&>* enumerate_neighbors_cyclic(
            ObjectContainer<sphere<double>>& cntnr, const typename mapped_type& sphere)
        {
            return util::make_enumerator(generator_type(
                boost::bind(&gen_take_neighbor_cyclic, _1, cntnr, sphere)));
        }

    private:
        Generators() {}

        inline static const result_type* gen_take_neighbor(
                generator_type::self& self,
                ObjectContainer<sphere<double>>& cntnr,
                const typename mapped_type& sphere)
        {
            collector col(self, cntnr.impl_.end());
            ::take_neighbor(cntnr.impl_, col, sphere);
            self.yield(NULL);
            self.exit();
            return NULL; // never get here
        }

        inline static const result_type* gen_take_neighbor_cyclic(
                generator_type::self& self,
                ObjectContainer<sphere<double>>& cntnr,
                const typename mapped_type& sphere)
        {
            collector col(self, cntnr.impl_.end());
            ::take_neighbor_cyclic(cntnr.impl_, col, sphere);
            self.yield(NULL);
            self.exit();
            return NULL; // never get here
        }
    };
#endif /* OBJECTMATRIX_USE_ITERATOR */

    class Builders
    {
    public:
        typedef std::vector<typename impl_type::iterator> sphere_ref_array_type;
        typedef std::vector<typename impl_type::key_type> key_array_type;
        typedef std::vector<double, util::pyarray_backed_allocator<double> >
                distance_array_type;
        //typedef boost::tuple<sphere_ref_array_type, distance_array_type>
//                result_type;
        typedef boost::tuple<key_array_type, distance_array_type>
                result_type;

        struct collector: public std::binary_function<
                typename impl_type::reference, typename position_type::value_type, void>
        {
            //// These 2 typedefs are not used.
            typedef typename impl_type::iterator first_argument_type;
            typedef typename position_type::value_type second_argument_type;
            typedef void result_type;
        public:
            //// Used in build_neighbors_array_cyclic.
            //// Difference between this collector and the 
            //// all_neighbors_collector is that this one doesn't do any 
            //// distance calculations, but just collects.
            inline collector(typename Builders::result_type& result)
                : //sa_(boost::get<0>(result)),
                ka_(boost::get<0>(result)),
                da_(boost::get<1>(result)) {}

            //// This operator is called from cyclic_neighbor_filter in 
            //// filters.hpp. Arguments are the i-th object in the matrix and 
            //// the distance d towards it (calculated in 
            //// cyclic_neighbor_filter).
            inline void operator()(typename impl_type::iterator i,
                    const typename position_type::value_type& d)
            {
                //sa_.push_back(i);
                ka_.push_back((*i).first);
                //// Distance already calculated, just store it.
                da_.push_back(d);
            }

        private:
            //sphere_ref_array_type& sa_;
            key_array_type& ka_;
            distance_array_type& da_;
        };

        struct all_neighbors_collector: public std::binary_function<
                typename impl_type::reference, typename position_type::value_type, void>
        {
            typedef typename impl_type::iterator first_argument_type;
            typedef typename position_type::value_type second_argument_type;
            typedef void result_type;
        public:
            //// Note: position is an argument of constructor.
            inline all_neighbors_collector(typename Builders::result_type& result,
                    const position_type& pos)
                : //sa_(boost::get<0>(result)),
                ka_(boost::get<0>(result)),
                      da_(boost::get<1>(result)),
                      pos_(pos) {}

            //// This operator is used from each_neighbor_loops in 
            //// object_container.hpp.
            inline void operator()(typename impl_type::iterator i)
            {
                //sa_.push_back(i);
                ka_.push_back((*i).first);
                mapped_type object = (*i).second;
                // Calculate and store distance from pos_ to object.
                da_.push_back(object.calculateDistanceToSelf(pos_));

            }

            //// This operator is used from each_neighbor_*cyclic*_loops in 
            //// object_container.hpp.
            inline void operator()(typename impl_type::iterator i,
                    const position_type& d)
            {
                //sa_.push_back(i);
                ka_.push_back((*i).first);
                //// Here finally the distance to the neighbor is calculated.
                //// The 'd' somehow takes care of cyclic boundaries I think.
                mapped_type object = (*i).second;
                da_.push_back(object.calculateDistanceToSelfWithOffset(pos_, d));
            }

        private:
            //sphere_ref_array_type& sa_;
            key_array_type& ka_;
            distance_array_type& da_;
            position_type pos_;
        };


    public:
        inline static void
        build_neighbors_array(result_type& retval,
                              ObjectContainer<mapped_type>& cntnr, const sphere_type& sphere)
        {
            //// Collect all objects from cntnr that lie within sphere.
            collector col(retval);
            take_neighbor(cntnr.impl_, col, sphere);
        }

        //// Called via some steps from  
        //// subSpaceSimulator.getNeighborsWithinRadiusNoSort().
        inline static void
        build_neighbors_array_cyclic(result_type& retval,
                ObjectContainer<mapped_type>& cntnr, const sphere_type& sphere)
        {
            //// Note: collector, not all_neighbors_collector.
            collector col(retval);
            //// This take_neighbor_cyclic is implemented in filters.hpp.
            take_neighbor_cyclic(cntnr.impl_, col, sphere);
        }

        inline static void
        build_all_neighbors_array(result_type& retval,
                ObjectContainer<mapped_type>& cntnr, const position_type& pos)
        {
            all_neighbors_collector col(retval, pos);
            cntnr.impl_.each_neighbor(cntnr.impl_.index(pos), col);
        }

        //// Called via some steps from subSpaceSimulator.getNeighbors().
        inline static void
        build_all_neighbors_array_cyclic(result_type& retval,
                ObjectContainer<mapped_type>& cntnr, const position_type& pos)
        {
            all_neighbors_collector col(retval, pos);
            //// cntnr refers to *this from all_neighbors_array_cyclic,
            //// so cntnr.impl_ is a implemented object_container from 
            //// object_container.hpp.
            cntnr.impl_.each_neighbor_cyclic(cntnr.impl_.index(pos), col);
        }

    private:
        Builders() {}
    };

public:
    ObjectContainer() {}

    ObjectContainer(length_type world_size, matrix_size_type size)
        : impl_(world_size, size) {}

    size_type size() const
    {
        return impl_.size();
    }

    size_type matrix_size() const
    {
        return impl_.matrix_size();
    }

    length_type world_size() const
    {
        return impl_.world_size();
    }

    length_type cell_size() const
    {
        return impl_.cell_size();
    }

#if OBJECTMATRIX_USE_ITERATOR
    //// Todo.
    //// When converting from spheres to cylinders, don't care about this part 
    //// because we're not using this iterneighbors.
    Enumerator<const Generators::result_type&>* iterneighbors(
            const sphere<double>& sphere)
    {
        return Generators::enumerate_neighbors(*this, sphere);
    }

    Enumerator<const Generators::result_type&>* iterneighbors_cyclic(
            const sphere<double>& sphere)
    {
        return Generators::enumerate_neighbors_cyclic(*this, sphere);
    }
#endif /* OBJECTMATRIX_USE_ITERATOR */

    boost::shared_ptr<typename Builders::result_type>
    //neighbors_array(const sphere& sphere)
    neighbors_array(const position_type& pos, const double radius)
    {
        typename Builders::distance_array_type::allocator_type alloc;

        boost::shared_ptr<typename Builders::result_type> retval(
            new typename Builders::result_type(
                typename boost::tuples::element<0, typename Builders::result_type>::type(),
                typename boost::tuples::element<1, typename Builders::result_type>::type(alloc)));
        Builders::build_neighbors_array(*retval, *this,
                                        sphere_type( pos, radius ) );

        // take over the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    //// Called via some steps from  
    //// subSpaceSimulator.getNeighborsWithinRadiusNoSort().
    //// This finds all neighbors within radius + distances towards them.  The 
    //// 'sphere_type' mentioned here is just the combination of the given 
    //// position and radius, and has nothing to do with the type of objects 
    //// that are stored in the matrix (which can also be cylinders for 
    //// example).
    boost::shared_ptr<typename Builders::result_type>
    //neighbors_array_cyclic(const sphere& sphere)
    neighbors_array_cyclic(const position_type& pos, const double radius)
    {
        typename Builders::distance_array_type::allocator_type alloc;

        boost::shared_ptr<typename Builders::result_type> retval(
            new typename Builders::result_type(
                typename boost::tuples::element<0, typename Builders::result_type>::type(),
                typename boost::tuples::element<1, typename Builders::result_type>::type(alloc)));
        Builders::build_neighbors_array_cyclic(*retval, *this,
                                               sphere_type( pos, radius ) );

        // take over the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    //// Note prescript all_.
    boost::shared_ptr<typename Builders::result_type>
    all_neighbors_array(const position_type& pos)
    {
        typename Builders::distance_array_type::allocator_type alloc;

        boost::shared_ptr<typename Builders::result_type> retval(
            new typename Builders::result_type(
                typename boost::tuples::element<0, typename Builders::result_type>::type(),
                typename boost::tuples::element<1, typename Builders::result_type>::type(alloc)));
        Builders::build_all_neighbors_array(*retval, *this, pos);

        // take over the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    //// Note prescript all_.
    //// This finds all neighbors + the distance towards them (compare 
    //// neighbors_array_cyclic) .
    //// Called via some steps from subSpaceSimulator.getNeighbors().
    boost::shared_ptr<typename Builders::result_type>
    all_neighbors_array_cyclic(const position_type& pos)
    {
        typename Builders::distance_array_type::allocator_type alloc;

        boost::shared_ptr<typename Builders::result_type> retval(
            new typename Builders::result_type(
                typename boost::tuples::element<0, typename Builders::result_type>::type(),
                typename boost::tuples::element<1, typename Builders::result_type>::type(alloc)));
        Builders::build_all_neighbors_array_cyclic(*retval, *this, pos);

        // take over the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    const bool contains( const key_type& k )
    {
        typename impl_type::iterator i(impl_.find(k));
        if (i == impl_.end())
        {
            return false;
        }
        return true;
    }

    void erase(const key_type& key)
    {
        impl_.erase(key);
    }

    operator impl_type&()
    {
        return impl_;
    }

    operator const impl_type&() const
    {
        return impl_;
    }

    //// Pythonify.
    inline static void __register_class()
    {
        using namespace boost::python;

        util::register_tuple_converter<boost::tuple<position_type,double> >();

#if OBJECTMATRIX_USE_ITERATOR
        util::register_tuple_converter<Generators::result_type>();

        util::register_enumerator<
                const Generators::result_type&,
                return_value_policy<copy_const_reference> >(
                "ObjectContainer_NeighborIterator");
#endif /* OBJECTMATRIX_USE_ITERATOR */

        util::register_multi_array_converter<
            typename boost::tuples::element<0, typename Builders::result_type>::type>();
        util::register_multi_array_converter<
            typename boost::tuples::element<1, typename Builders::result_type>::type>();
        // the following conversion is the same as the previous
        // util::register_multi_array_converter<boost::tuples::element<2, Builders::result_type>::type>();

        util::register_tuple_converter<typename Builders::result_type>();
    }

    inline static void __register_object_container()
    {
        using namespace boost::python;

        typedef ObjectContainer< mapped_type > ObjectContainerType;
        class_<ObjectContainerType>("ObjectContainer")
        .def(init<typename ObjectContainerType::length_type, typename ObjectContainerType::matrix_size_type>())
#if OBJECTMATRIX_USE_ITERATOR
        .def("iterneighbors", &ObjectContainerType::iterneighbors,
                return_value_policy<manage_new_object>())
        .def("iterneighbors_cyclic", &ObjectContainerType::iterneighbors_cyclic,
                return_value_policy<manage_new_object>())
#endif  /* OBJECTMATRIX_USE_ITERATOR */
        
        .add_property("cell_size", &ObjectContainerType::cell_size)
        .add_property("world_size", &ObjectContainerType::world_size)
        .add_property("matrix_size", &ObjectContainerType::matrix_size)
        .def("neighbors_array", &ObjectContainerType::neighbors_array)
        .def("neighbors_array_cyclic", &ObjectContainerType::neighbors_array_cyclic)
        .def("all_neighbors_array", &ObjectContainerType::all_neighbors_array)
        .def("all_neighbors_array_cyclic", &ObjectContainerType::all_neighbors_array_cyclic)
        .def("size", &ObjectContainerType::size)
        .def("contains", &ObjectContainerType::contains)
        .def("erase", &ObjectContainerType::erase);
    }

protected:
    //// Implementation type of an object container, once all templates are 
    //// specified.
    impl_type impl_;
};


class SphereContainer: public ObjectContainer< sphere<double> >
{
public:
    SphereContainer() {}

    SphereContainer(length_type world_size, matrix_size_type size)
        : ObjectContainer<mapped_type>(world_size, size) {}

    const boost::tuple<position_type,double> get( const key_type& k )
    {
        impl_type::iterator i(impl_.find(k));
        if (i == impl_.end())
        {
            throw std::runtime_error( "key not found." );
        }

        return boost::make_tuple( (*i).second.position, 
                                  (*i).second.radius );
    }


    void insert( const key_type& key, const position_type& p, const double r )
    {
        //// Called from cObjectMatrix.py.
        //// Calls insert() in object_container.hpp, after first making an 
        //// instance of the value_type, which has arguments:
        //// + key type: a python object (single/pair).
        //// + mapped type: a *sphere* at position p with radius r.
        impl_.insert(impl_type::value_type(key, mapped_type(p,r)));
    }


    void update( const key_type& key, const position_type& p, const double r )
    {
        impl_.insert(impl_type::value_type(key, mapped_type(p,r)));
    }


    inline static void __register_class()
    {
        peer::ObjectContainer<mapped_type>::__register_object_container();

        using namespace boost::python;
        class_<SphereContainer, bases<ObjectContainer<mapped_type> > >("SphereContainer")
            .def(init<SphereContainer::length_type, SphereContainer::matrix_size_type>())
            .def("insert", &SphereContainer::insert)
            .def("update", &SphereContainer::update)
            .def("get", &SphereContainer::get);
    
    }
};


class CylinderContainer: public ObjectContainer< cylinder<double> >
{
public:
    CylinderContainer() {}

    CylinderContainer(length_type world_size, matrix_size_type size)
        : ObjectContainer< mapped_type >(world_size, size) {}

    void insert( const key_type& key, const position_type& position, const position_type& orientation, const double radius, const double halfLength )
    {
        impl_.insert(impl_type::value_type(key, mapped_type(position, orientation, radius, halfLength)));
    }

    void update( const key_type& key, const position_type& position, const position_type& orientation, const double radius, const double halfLength )
    {
        impl_.insert(impl_type::value_type(key, mapped_type(position, orientation, radius, halfLength)));
    }

    inline static void __register_class()
    {
        peer::ObjectContainer<mapped_type>::__register_object_container();

        using namespace boost::python;
        class_<CylinderContainer, bases<ObjectContainer<mapped_type> > >("CylinderContainer")
            .def(init<CylinderContainer::length_type, CylinderContainer::matrix_size_type>())
            .def("insert", &CylinderContainer::insert)
            .def("update", &CylinderContainer::update);
    }
};


}// namespace peer

#endif /* OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP */
