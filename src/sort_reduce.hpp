/***
 *  $Id: sort_reduce.hpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: sort_reduce.hpp
 *  Created: Apr 18, 2008
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

#ifndef SORT_REDUCE_HPP
#define SORT_REDUCE_HPP

#include <string.h>
#include <mpi.h>


template <typename T, template <typename U> class Op>
void sort_reduce(void* in, void* inout, int* len, MPI_Datatype*) {
    int n = *len;
    T* buf = new T[n];

    T* pos1 = static_cast<T*>(in);
    T* pos2 = static_cast<T*>(inout);

    Op<T> comp;

    for (int i = 0; i < n; ++i) {
	if (comp(*pos1, *pos2) == true) buf[i] = *(pos1++);
	else if (comp(*pos2, *pos1) == true) buf[i] = *(pos2++);
	else {
	    buf[i] = *pos1;
	    ++pos1;
	    ++pos2;
	}
    }

    memcpy(inout, buf, n * sizeof(T));

    delete[] buf;
} // sort_reduce

#endif // SORT_REDUCE_HPP
