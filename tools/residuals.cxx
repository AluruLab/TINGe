/***
 *  $Id: residuals.cxx 597 2010-06-23 01:06:56Z zola $
 **
 *  File: residuals.cxx
 *  Created: Feb 12, 2009
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

#include <mpix/MPI_env.hpp>
#include <jaz/stat_util.hpp>
#include <jaz/sys_tools.hpp>
#include "AppConfig.hpp"
#include "InputData.hpp"
#include "read_exp.hpp"
#include "write_block.hpp"


int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(0);

    MPI_Init(&argc, &argv);

    // create 1D topology
    int wsize = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);

    int period = 1;

    MPI_Comm mpi_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 1, &wsize, &period, 1, &mpi_comm);

    // create MPI environment
    mpix::MPI_env mpi_env(mpi_comm);
    mpi_env.root(0);

    mpi_iom.mpi_env = &mpi_env;

    if (jaz::iss_pubsetbuf_test() == false) {
	std::cout << master << oops
		  << "your compiler does not support pubsetbuf(...) method!\n"
		  << "Change compiler or contact authors for advice."
		  << std::endl;

	return MPI_Finalize();
    }

    if (argc != 3) {
	std::cout << master << "usage: " << argv[0] << " infile outfile\n";
	return MPI_Finalize();
    }

    AppConfig app_conf;
    app_conf.input = argv[1];

    // read input data
    InputData data;
    bool res = rw_read_exp(mpi_env, app_conf, data);

    if (res == false) return MPI_Finalize();

    if (mpi_env.rank() == 0) {
	std::ostringstream os;
	for (unsigned int i = 0; i < data.g_ncol; ++i) {
	    os.str("");
	    os << "V" << (i + 1);
	    data.exp_table.col_name(i, os.str());
	}
    }

    data.exp_table.row_names(data.probes.begin() + data.first_row,
			     data.probes.begin() + data.last_row);

    std::cout << master << "* calculating residuals..." << std::endl;
    double t0 = MPI_Wtime();

    // here we go
    for (unsigned int i = 0; i < data.num_rows; ++i) {
	double m = jaz::mean<double>(data.exp_table.row_begin(i),
				     data.exp_table.row_end(i));
	double sse = 0;

	for (unsigned int j = 0; j < data.g_ncol; ++j) {
	    double x = data.exp_table[i][j] - m;
	    data.exp_table[i][j] = x;
	    sse += x * x;
	}

	sse = sse * (1 - (1.0 / data.g_ncol));

	for (unsigned int j = 0; j < data.g_ncol; ++j) {
	    double x = data.exp_table[i][j];
	    data.exp_table[i][j] = x * sqrt(data.g_ncol / (sse - (x * x)));
	}
    } // for i

    double t1 = MPI_Wtime();
    std::cout << master << "calculations done: " << timer(t1 - t0) << std::endl;

    // write deleted residuals
    std::cout << master << "* writing table " << argv[2] << " ..." << std::endl;
    t0 = MPI_Wtime();

    std::string s;
    std::stringstream ss;

    ss << data.exp_table;
    if (mpi_env.rank() != 0) std::getline(ss, s);

    write_block(mpi_env, argv[2], ss.str().c_str(), ss.str().size());

    t1 = MPI_Wtime();
    std::cout << master << "writing done: " << timer(t1 - t0) << std::endl;

    return MPI_Finalize();
} // main
