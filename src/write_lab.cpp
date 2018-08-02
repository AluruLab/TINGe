/***
 *  $Id: write_lab.cpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: write_lab.cpp
 *  Created: Apr 20, 2008
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

#include "write_lab.hpp"


bool rw_write_lab(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
		  const InputData& data,
		  const jaz::plain_array<double>& adj_matrix,
		  const std::string& name) {
    std::cout << master << "* writing MCL label file " << name << " ..."
	      << std::endl;
    TIMER_START;

    std::ostringstream ss;

    // write adj matrix
    double val;
    for (unsigned int i = 0; i < data.num_rows; ++i) {
	for (unsigned int j = i + data.first_row + 1; j < data.g_nrow; ++j) {
	    val = adj_matrix[i * data.g_nrow + j];
	    if (0.0 < val) {
		ss << data.probes[data.first_row + i] << ' ' << data.probes[j]
		   << ' ' << val << '\n';
	    }
	}
    }

    if (app_conf.mem_report == true) MEM_REPORT;

    unsigned int size = ss.str().size();
    write_block(mpi_env, name.c_str(), ss.str().c_str(), size);

    std::cout << master << "writing done: " << timer(TIMER_GET) << std::endl;

    return true;
} // rw_write_lab

