/***
 *  $Id: transform_data.cpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: rank_transform.cpp
 *  Created: Mar 29, 2008
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

#include "transform_data.hpp"


void sample_data(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
		 InputData& data) {

    std::cout << master << "* sampling data..." << std::endl;

    if (data.g_ncol < app_conf.sample_size) {
	std::cout << master << warning
		  << "sample size too large, sampling ignored" << std::endl;
	return;
    }

    unsigned int m = app_conf.sample_size;

    // sample data
    jaz::plain_array<unsigned int> id(data.g_ncol);
    for (unsigned int i = 0; i < data.g_ncol; i++) id[i] = i;

    jaz::plain_array<unsigned int> sample(app_conf.sample_size);

    if (mpi_env.am_I_root() == true) {
	jaz::random_sample_n(id.begin(), id.end(), sample.begin(), m);

	std::string name = app_conf.output;
	std::string::size_type pos = name.rfind('.');

	if (pos != std::string::npos) {
	    name = std::string(name.begin(), name.begin() + pos);
	}

	name += ".sample";

	std::ofstream f(name.c_str());
	std::copy(sample.begin(), sample.end(),
		  std::ostream_iterator<unsigned int>(f, "\n"));
	f.close();

	std::cout << "column ids written to " << name << std::endl;
    } // if root

    MPI_Bcast(sample.begin(), m, MPI_UNSIGNED, mpi_env.root(), mpi_env.comm());

    // extract sample
    jaz::numeric_table<double> exp_table(data.num_rows, m);

    if (app_conf.mem_report == true) MEM_REPORT;

    for (unsigned int i = 0; i < data.num_rows; ++i) {
	double* tab = data.exp_table[i];
	double* dest = exp_table[i];
	for (unsigned int j = 0; j < m; ++j) dest[j] = tab[sample[j]];
    }

    data.exp_table = exp_table;
    data.g_ncol = m;
} // sample_data

void rank_transform(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
		    InputData& data) {

    std::cout << master << "* rank transforming expression data" << std::endl;

    unsigned int n = data.num_rows;
    unsigned int ncol = data.g_ncol;

    data.rank_array.resize(n * ncol);

    for (unsigned int i = 0; i < n; ++i) {
	jaz::rank_transform(data.exp_table.row_begin(i), data.exp_table.row_end(i),
			    data.rank_array.begin() + i * ncol);
    }
} // rank_transform
