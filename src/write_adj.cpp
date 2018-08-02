/***
 *  $Id: write_adj.cpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: write_adj.cpp
 *  Created: Apr 08, 2008
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

#include "write_adj.hpp"


bool rw_write_adj(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
		  const InputData& data,
		  const jaz::plain_array<double>& adj_matrix,
		  const std::string& name) {
    std::cout << master << "* writing network " << name << " ..." << std::endl;
    TIMER_START;

    // get network properties
    unsigned int max_edg = ((data.g_nrow * data.g_nrow) - data.g_nrow);
    unsigned int tot_edg = get_tot_edges(mpi_env, data, adj_matrix);

    double fdr = (data.fdr < 0.0) ?
	(app_conf.mi_pval * max_edg) / (2.0 * tot_edg) : data.fdr;

    double d = (2.0 * tot_edg) / max_edg;

    std::ostringstream ss;

    // write ADJ header
    if (mpi_env.am_I_root() == true) {
	const int N = 32;

	ss.flags(std::ios_base::left);
	ss << std::setfill('.');

	ss << std::setw(N) << "> TINGe version: " << ' '
	   << SHORT_NAME << " " << MAJOR_VERSION << "." << MINOR_VERSION << '\n';

	std::time_t t = std::time(&t);
	ss << std::setw(N) << "> Date: " << ' ' << std::ctime(&t);

	ss << std::setw(N) << "> Cmd line: ";
	for (int i = 0; i < app_conf.argc; ++i) {
	    ss << ' ' << app_conf.argv[i];
	}
	ss << '\n';

	ss << std::setw(N) << "> CPUs used: " << ' ' << mpi_env.size() << '\n';

	ss << std::setw(N) << "> Input size: " << ' ' << data.g_nrow;
	ss << "x" << data.g_ncol << '\n';

	ss << std::setw(N) << "> Execution time: " << ' ' << data.etime;
	ss << " [s]\n";

	app_conf.print_config(ss);

	ss << std::setw(N) << "> Random MI, mean: " << ' ' << data.nd_mean << '\n';
	ss << std::setw(N) << "> Random MI, stdev: " << ' ' << data.nd_sd << '\n';

	ss << std::setw(N) << "> MI threshold used: " << ' ' << data.mi_thr_used << '\n';

	ss << std::setw(N) << "> Final MI, mean: " << ' ' << data.mi_mean << '\n';
	ss << std::setw(N) << "> Final MI, stdev: " << ' ' << data.mi_sd << '\n';

	ss << std::setw(N) << "> Network dimensions: " << ' '
	   << data.g_nrow << "x" << tot_edg << '\n';

	ss << std::setw(N) << "> Network density: " << ' ' << d << '\n';

	ss << std::setw(N) << "> Estimated FDR: " << ' ' << fdr << '\n';

	if (data.boot > 0) {
	    ss << std::setw(N) << "> Bootstrapped network: " << ' '
	       << data.boot << '\n';
	}
    } // if mpi_env

    // write adj matrix
    double val;
    for (unsigned int i = 0; i < data.num_rows; ++i) {
	unsigned int j = i + data.first_row + 1;

	for (; j < data.g_nrow; ++j) {
	    val = adj_matrix[i * data.g_nrow + j];
	    if (0.0 < val) break;
	}

	if (j < data.g_nrow) {
	    ss << data.probes[data.first_row + i];
	    for (j = i + data.first_row + 1; j < data.g_nrow; ++j) {
		val = adj_matrix[i * data.g_nrow + j];
		if (0.0 < val) ss << '\t' << data.probes[j] << '\t' << val;
	    }
	    ss << '\n';
	}
    }

    if (app_conf.mem_report == true) MEM_REPORT;

    unsigned int size = ss.str().size();
    write_block(mpi_env, name.c_str(), ss.str().c_str(), size);

    std::cout << master << "writing done: " << timer(TIMER_GET) << std::endl;

    return true;
} // rw_write_adj
