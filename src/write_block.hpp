/***
 *  $Id: write_block.hpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: write_block.hpp
 *  Created: Apr 21, 2008
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

#ifndef WRITE_BLOCK_HPP
#define WRITE_BLOCK_HPP

#include <mpix/MPI_env.hpp>
#include <jaz/sys_tools.hpp>


inline bool write_block(const mpix::MPI_env& mpi_env,
			const char* name, const char* buf, unsigned int n) {
    if (jaz::file_size(name) != 0) {
	MPI_File_delete(const_cast<char*>(name), MPI_INFO_NULL);
    }

    MPI_File fh;
    MPI_Status stat;

    MPI_File_open(mpi_env.comm(), const_cast<char*>(name),
		  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    unsigned int size = n;
    unsigned int offset = 0;

    MPI_Scan(&size, &offset, 1, MPI_UNSIGNED, MPI_SUM, mpi_env.comm());
    offset -= size;

    char type[] = "native";

    MPI_File_set_view(fh, offset, MPI_CHAR, MPI_CHAR, type, MPI_INFO_NULL);
    MPI_File_write_all(fh, const_cast<char*>(buf), size, MPI_CHAR, &stat);
    MPI_File_close(&fh);

    return true;
} // write_block

#endif // WRITE_BLOCK_HPP
