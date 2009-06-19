#ifndef OBJECT_CONTAINER_HPP
#define OBJECT_CONTAINER_HPP

#include <cstddef>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/list/at.hpp>
#include "array_helper.hpp"
#include "vmap.hpp"

#include "utils.hpp"


//// This template is called from ../boost.python/peer/ObjectContainer.hpp 
//// with the arguments: <double, key_type, mapped_type, get_mapper_mf>, and 
//// is called impl_type then.
//// MFget_mapper is from now on the type for a mapper, like for example 
//// std::map.
////
//// This class seems to be competely generic. So for example nowhere it is 
//// assumed the mapped_type is a necessarily a sphere.
template<typename T_, typename Tkey_, typename Tmapped_type_,
        template<typename, typename> class MFget_mapper_ =
            get_default_impl::std::template map>
class object_container
{
public:
    //// length_type
    //// For example micrometers.
    typedef T_ length_type;

    //// key_type
    //// For example a python object.
    typedef Tkey_ key_type;

    typedef Tmapped_type_ mapped_type;

    //// value_type
    //// Pair of a key and a sphere.
    typedef std::pair<const key_type, mapped_type> value_type;

    //// position_type
    typedef position<length_type> position_type;

    //// cell_type
    ////   + key: python object (in this case)
    ////   + mapped_type: sphere (in this case)
    ////   + MVget_mapper: mapper (see ObjectContainer.hpp)
    typedef vmap<Tkey_, mapped_type, MFget_mapper_> cell_type;

    //// matrix_type
    //// Matrix is a 3 dimensional array of type cell_type.
    typedef boost::multi_array<cell_type, 3> matrix_type;

    //// size_type ...
    typedef typename cell_type::size_type size_type;

    //// cell_index_type ...
    typedef boost::array<typename matrix_type::size_type, 3>
            cell_index_type;

    //// cell_offset_type ...
    typedef boost::array<typename matrix_type::difference_type, 3>
            cell_offset_type;

    typedef typename cell_type::reference reference;

    typedef typename cell_type::const_reference const_reference;

    typedef typename MFget_mapper_<key_type, cell_type*>::type
            id_to_cell_mapper_type;

    typedef id_to_cell_mapper_type key_to_cell_mapper_type;




    //// Iterator base.
    template<typename Thost_, typename Treftype_, typename Tcreftype_,
            typename Tciter_>
    class iterator_base
    {
    public:
        typedef Treftype_ reference;
        typedef Tcreftype_ const_reference;

    public:
        iterator_base(cell_type* cell_p, const Tciter_& cell_iter_end,
                const Tciter_& cell_iter)
                : cell_p_(cell_p), cell_iter_end_(cell_iter_end),
                  cell_iter_(cell_iter) {}

        bool operator==(const Thost_& rhs)
        {
            return cell_p_ == rhs.cell_p_ &&
                cell_iter_ == rhs.cell_iter_;
        }

        bool operator!=(const Thost_& rhs)
        {
            return !operator==(rhs);
        }

        const_reference operator*() const
        {
            return *cell_iter_;
        }

        reference operator*()
        {
            return *cell_iter_;
        }

        Thost_& operator=(const iterator_base& rhs)
        {
            cell_p_ = rhs.cell_p_;
            cell_iter_end_ = rhs.cell_iter_end_;
            cell_iter_ = rhs.cell_iter_;
            return *reinterpret_cast<Thost_*>(this);
        }

    public:
        cell_type* cell_p_;
        mutable Tciter_ cell_iter_end_;
        Tciter_ cell_iter_;
    };



    //// Iterator.
    class iterator: public iterator_base<iterator, reference, const_reference,
            typename cell_type::iterator>
    {
    public:
        typedef iterator_base<iterator, reference, const_reference,
                typename cell_type::iterator> base_type;

    public:
        iterator(cell_type* cell_p,
            const typename cell_type::iterator& cell_iter_end,
            const typename cell_type::iterator& cell_iter)
                : base_type(cell_p, cell_iter_end, cell_iter) {}
    };


    //// Const iterator base.
    class const_iterator: public iterator_base<const_iterator,
            const_reference, const_reference,
            typename cell_type::const_iterator>
    {
    public:
        typedef iterator_base<const_iterator, const_reference, const_reference,
                typename cell_type::const_iterator> base_type;

    public:
        const_iterator(cell_type* cell_p,
            const typename cell_type::const_iterator& cell_iter_end,
            const typename cell_type::const_iterator& cell_iter)
                : base_type(cell_p, cell_iter_end, cell_iter) {}

        const_iterator(const iterator& that)
                : base_type(that.cell_p_, that.cell_iter_end_, that.cell_iter_) {}
    };

public:
    //// Constructor.
    //// Two paramaters: worldSize and matrixSize.
    object_container(length_type world_size = 1.0,
            typename matrix_type::size_type size = 1)
        : world_size_(world_size),
          cell_size_(world_size / size),
          matrix_(boost::extents[size][size][size]),
          size_(0)
    {
    }


    //// Called from insert().
    inline cell_index_type index(const position_type& pos,
            double t = 1e-10) const
    {
        return array_gen<typename matrix_type::size_type>(
            static_cast<typename matrix_type::size_type>(
                trunc( pos.x() / cell_size_ ) ) % matrix_.shape()[0],
            static_cast<typename matrix_type::size_type>(
                trunc( pos.y() / cell_size_ ) ) % matrix_.shape()[1],
            static_cast<typename matrix_type::size_type>(
                trunc( pos.z() / cell_size_ ) ) % matrix_.shape()[2] );
    }


    //// Called from other file?
    inline cell_offset_type offset(const position_type& pos,
            double t = 1e-10) const
    {
        return array_gen<typename matrix_type::difference_type>(
            trunc( pos.x() / cell_size ),
            trunc( pos.y() / cell_size ),
            trunc( pos.z() / cell_size ) );
    }


    //// Called by each_neighbor_loops().
    inline bool offset_index(
            cell_index_type& i,
            const cell_offset_type& o) const
    {
        if ((o[0] < 0 && static_cast<size_type>(-o[0]) > i[0])
                || (matrix_.shape()[0] - o[0] <= i[0])
                || (o[1] < 0 && static_cast<size_type>(-o[1]) > i[1])
                || (matrix_.shape()[1] - o[1] <= i[1])
                || (o[2] < 0 && static_cast<size_type>(-o[2]) > i[2])
                || (matrix_.shape()[2] - o[2] <= i[2]))
        {
            return false;
        }
        i[0] += o[0];
        i[1] += o[1];
        i[2] += o[2];
        return true;
    }


    inline position_type offset_index_cyclic(cell_index_type& i,
                                             const cell_offset_type& o)
    {
        position_type retval;

        if (o[0] < 0 &&
            static_cast<typename matrix_type::size_type>(-o[0]) > i[0])
        {
            typename matrix_type::size_type t(
                (i[0] + matrix_.shape()[0] - (-o[0] % matrix_.shape()[0])) %
                matrix_.shape()[0]);
            retval[0] 
                = (o[0] - 
                   static_cast<typename matrix_type::difference_type>
                   (t - i[0])) * cell_size_;
            i[0] = t;
        }
        else if (matrix_.shape()[0] - o[0] <= i[0])
        {
            typename matrix_type::size_type t(
                    (i[0] + (o[0] % matrix_.shape()[0])) % matrix_.shape()[0]);
            retval[0] 
                = (o[0] - 
                   static_cast<typename matrix_type::difference_type>
                   (t - i[0])) * cell_size_;
            i[0] = t;
        }
        else
        {
            i[0] += o[0];
        }

        if (o[1] < 0 &&
                static_cast<typename matrix_type::size_type>(-o[1]) > i[1])
        {
            typename matrix_type::size_type t(
                    (i[1] + matrix_.shape()[1] - (-o[1] % matrix_.shape()[1])) %
                        matrix_.shape()[1]);
            retval[1] = (o[1] - static_cast<typename matrix_type::difference_type>(t - i[1])) * cell_size_;
            i[1] = t;
        }
        else if (matrix_.shape()[1] - o[1] <= i[1])
        {
            typename matrix_type::size_type t(
                    (i[1] + (o[1] % matrix_.shape()[1])) % matrix_.shape()[1]);
            retval[1] = (o[1] - static_cast<typename matrix_type::difference_type>(t - i[1])) * cell_size_;
            i[1] = t;
        }
        else
        {
            i[1] += o[1];
        }

        if (o[2] < 0 &&
                static_cast<typename matrix_type::size_type>(-o[2]) > i[2])
        {
            typename matrix_type::size_type t(
                    (i[2] + matrix_.shape()[2] - (-o[2] % matrix_.shape()[2])) %
                        matrix_.shape()[2]);
            retval[2] = (o[2] - static_cast<typename matrix_type::difference_type>(t - i[2])) * cell_size_;
            i[2] = t;
        }
        else if (matrix_.shape()[2] - o[2] <= i[2])
        {
            typename matrix_type::size_type t(
                    (i[2] + (o[2] % matrix_.shape()[2])) % matrix_.shape()[2]);
            retval[2] = (o[2] - static_cast<typename matrix_type::difference_type>(t - i[2])) * cell_size_;
            i[2] = t;
        }
        else
        {
            i[2] += o[2];
        }

        return retval;
    }


    inline const cell_type& cell(const cell_index_type& i) const
    {
        return matrix_[i[0]][i[1]][i[2]];
    }


    //// A cell is an element of a matrix, indexed by a cell_index_type.
    inline cell_type& cell(const cell_index_type& i)
    {
        return matrix_[i[0]][i[1]][i[2]];
    }


    inline length_type world_size() const
    {
        return world_size_;
    }


    inline length_type cell_size() const
    {
        return cell_size_;
    }


    inline typename matrix_type::size_type matrix_size() const
    {
        return matrix_.shape()[0];
    }


    inline size_type size() const
    {
        return size_;
    }


    //// Called from ObjectContainer.hpp.
    //// The actual insert.
    inline std::pair<iterator, bool> insert(const value_type& v)
    {
        //// Construct a cell c.
        //// typedef vmap<Tkey_, mapped_type, MFget_mapper_> cell_type;
        //// Todo. I still don't understand this index business. 
        cell_type& c(cell(index(v.second.position)));
        typename key_to_cell_mapper_type::iterator kci(rmap_.find(v.first));
        if (rmap_.end() != kci)
        {
            if (&c != (*kci).second)
            {
                (*kci).second->erase(v.first);
                rmap_.erase(v.first);
                rmap_.insert(std::make_pair(v.first, &c));
            }
            std::pair<typename cell_type::iterator, bool> ir(c.insert(v));
            return std::make_pair(iterator(&c, c.end(), ir.first), false);
        }
        std::pair<typename cell_type::iterator, bool> ir(c.insert(v));
        rmap_.insert(std::make_pair(v.first, &c));
        ++size_;
        return std::make_pair(iterator(&c, c.end(), ir.first), ir.second);
    }


    //// Called by insert().
    inline bool erase(const key_type& k)
    {
        typename key_to_cell_mapper_type::iterator p(rmap_.find(k));
        if (rmap_.end() == p)
        {
            return false;
        }
        (*p).second->erase(k);
        rmap_.erase(p);
        --size_;
        return true;
    }


    //// Called by each_neighbor_loops() each_neighbour_cyclic_loop()
    inline iterator begin()
    {
        cell_type* cell_p(matrix_.origin());
        cell_type* cell_pe(matrix_.origin() + matrix_.num_elements());

        while (cell_p->begin() == cell_p->end() && cell_p < cell_pe)
        {
            ++cell_p;
        }

        return iterator(cell_p, cell_p->begin(), cell_p->begin());
    }


    //// Called many times.
    inline iterator end()
    {
        cell_type* last(matrix_.origin() + matrix_.num_elements() - 1);
        return iterator(last, last->end(), last->end());
    }


    //// Called by insert()/erase().
    inline iterator find(const key_type& k)
    {
        typename key_to_cell_mapper_type::iterator p(rmap_.find(k));
        if (rmap_.end() == p)
        {
            return end();
        }
        typename cell_type::iterator i((*p).second->find(k));
        if ((*p).second->end() == i)
        {
            return end();
        }
        return iterator((*p).second, (*p).second->end(), i);
    }


    //// Called by insert()/erase().
    inline const_iterator find(const key_type& k) const
    {
        typename key_to_cell_mapper_type::const_iterator p(rmap_.find(k));
        if (rmap_.end() == p)
        {
            return end();
        }
        typename cell_type::const_iterator i((*p).second->find(k));
        if ((*p).second->end() == i)
        {
            return end();
        }
        return const_iterator((*p).second, (*p).second->end(), i);
    }


    template<typename Tcollect_>
    inline void each_neighbor(const cell_index_type& idx, Tcollect_& collector)
    {
        cell_offset_type _off;
        each_neighbor_loops<Tcollect_&>(3, _off, idx, collector);
    }


    template<typename Tcollect_>
    inline void each_neighbor(const cell_index_type& idx,
                              const Tcollect_& collector)
    {
        cell_offset_type _off;
        each_neighbor_loops<const Tcollect_&>(3, _off, idx, collector);
    }


    //// Called via some steps from subSpaceSimulator.getNeighbors().
    //// But also used for subSpaceSimulator.getNeighborsWithinRadiusNoSort(), 
    //// but then called via take_neighbor_cyclic in filters.hpp.
    template<typename Tcollect_>
    inline void each_neighbor_cyclic(const cell_index_type& idx,
            Tcollect_& collector)
    {
        cell_offset_type _off;
        each_neighbor_cyclic_loops<Tcollect_&>(3, _off, idx, collector);
    }


    template<typename Tcollect_>
    inline void each_neighbor_cyclic(const cell_index_type& idx,
            const Tcollect_& collector)
    {
        cell_offset_type _off;
        each_neighbor_cyclic_loops<const Tcollect_&>(3, _off, idx, collector);
    }

private:

    template<typename Tcollect_>
    inline void each_neighbor_loops(const std::size_t depth,
            cell_offset_type& off, const cell_index_type& idx,
            Tcollect_ collector)
    {
        if (depth > 0)
        {
            for (off[depth - 1] = -1; off[depth - 1] <= 1; ++off[depth - 1])
            {
                each_neighbor_loops(depth - 1, off, idx, collector);
            }
        }
        else
        {
            cell_index_type _idx(idx);
            if (!offset_index(_idx, off)) {
                return;
            }
            cell_type& c(cell(_idx));
            for (typename cell_type::iterator i(c.begin()); i != c.end(); ++i) 
            {
                collector(iterator(&c, c.end(), i));
            }
        }
    }


    //// Called via some steps from subSpaceSimulator.getNeighbors().
    //// and subSpaceSimulator.getNeighborsWithinRadiusNoSort().
    template<typename Tcollect_>
    inline void each_neighbor_cyclic_loops(const std::size_t depth,
            cell_offset_type& off, const cell_index_type& idx,
            Tcollect_ collector)
    {
        if (depth > 0)
        {
            for (off[depth - 1] = -1; off[depth - 1] <= 1; ++off[depth - 1])
            {
                each_neighbor_cyclic_loops(depth - 1, off, idx, collector);
            }
        }
        else
        {
            cell_index_type _idx(idx);
            const position_type pos_off(offset_index_cyclic(_idx, off));
            cell_type& c(cell(_idx));
            for (typename cell_type::iterator i(c.begin()); i != c.end(); ++i) 
            {
                //// Here finally the distance to the neighbor is calculated.
                //// collector is all_neighbors_collector in 
                //// ObjectContainer.hpp for getNeighbors(), but 
                //// cyclic_neighbor_filter in filters.hpp for 
                //// getNeighborsWithinRadiusNoSort().
                collector(iterator(&c, c.end(), i), pos_off);
            }
        }
    }

private:
    const length_type world_size_;
    const length_type cell_size_;
    //// An object container has a matrix.
    matrix_type matrix_;
    id_to_cell_mapper_type rmap_;
    size_type size_;
};

//// Operator +=.
template<typename T_, typename Tkey_, typename mapped_type,
        template<typename, typename> class MFget_mapper_>
static inline typename object_container<T_, Tkey_, mapped_type, MFget_mapper_>::cell_index_type&
operator+=(
       typename object_container<T_,
                Tkey_, mapped_type, MFget_mapper_>::cell_index_type& lhs,
       const typename object_container<T_,
                Tkey_, mapped_type, MFget_mapper_>::cell_offset_type& rhs)
{
    rhs[0] += lhs[0];
    rhs[1] += lhs[1];
    rhs[2] += lhs[2];
    return rhs;
}


#endif /* OBJECT_CONTAINER_HPP */
