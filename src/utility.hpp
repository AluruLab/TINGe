/***
 *  $Id: utility.hpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: utility.hpp
 *  Created: Mar 28, 2008
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *
 *  Copyright 2007-2010 Jaroslaw Zola
 *
 *  This file is part of TINGe.
 *
 *  TINGe is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  TINGe is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with TINGe. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <jaz/sys_tools.hpp>
#include <utility>
#include <mpi.h>


#define TIMER_START double ts__ = MPI_Wtime()
#define TIMER_GET (MPI_Wtime() - ts__)


#define MEM_REPORT std::cout << master << "___memory usage: " \
    << ((1.0 * jaz::mem_usage()) / (1024 * 1024)) << " MB"    \
    << std::endl;


inline int row2cpu(int rank, int size, unsigned int nrow,
		   unsigned int& frow, unsigned int& lrow) {
    unsigned int n = nrow / size;
    unsigned int r = nrow % size;

    if (static_cast<unsigned int>(rank) < r) ++n;
    frow = rank * n;

    if (static_cast<unsigned int>(rank) >= r) frow += r;
    lrow = frow + n;

    return n;
} // row2cpu


// it is minimalized fake forward iterator to iterate over columns
template <typename T> class jump_iterator {
public:
    explicit jump_iterator(const T* from, unsigned int jump)
	: from_(from), jump_(jump) { }

    jump_iterator& operator++() {
	from_ += jump_;
	return *this;
    } // operator++

    const T& operator*() { return *from_; }

    bool operator==(const jump_iterator& iter) const {
	return (from_ == iter.from_);
    } // operator ==

    bool operator!=(const jump_iterator& iter) const {
	return (from_ != iter.from_);
    } // operator!=

private:
    const T* from_;
    unsigned int jump_;

}; // class jump_iterator


template <typename Iter, typename T>
std::pair<Iter, Iter> min_max_val(Iter first, Iter last, T val) {
    std::pair<Iter, Iter> res(last, last);

    while (first != last) {
	if (val < *first) {
	    res.first = first;
	    res.second = first;
	    break;
	}
	++first;
    }

    for (; first != last; ++first) {
	if (val < *first) {
	    if (*first < *res.first) res.first = first;
	    else if (*res.second < *first) res.second = first;
	}
    }

    return res;
} // min_max_val

#endif // UTILITY_HPP
