/***
 *  $Id: main.cxx 597 2010-06-23 01:06:56Z zola $
 **
 *  File: main.cxx
 *  Created: Nov 11, 2007
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
#include "AppConfig.hpp"
#include "InputData.hpp"
#include "adj_build.hpp"
#include "read_adj.hpp"
#include "read_exp.hpp"
#include "read_tfs.hpp"
#include "transform_data.hpp"
#include "write_adj.hpp"
#include "write_lab.hpp"


void check_tfs(InputData& data) {
    unsigned int atf = 0;
    data.isTF.resize(data.g_nrow, false);
    for (unsigned int i = 0; i < data.g_nrow; ++i) {
	std::string name = data.probes[i];
	if (data.tfs.find(name) != data.tfs.end()) {
	    data.isTF[i] = true;
	    atf++;
	}
    }
    if (atf > 0) {
	std::cout << master << "number of active TFs: " << atf << std::endl;
	data.has_TF = true;
    }
} // check_tfs


void check_test_size(AppConfig& app_conf, InputData& data) {
    // check if permutation size if sufficient (not too small, not too large)
    unsigned int n = data.g_nrow;

    unsigned int null_size
	= static_cast<unsigned int>(app_conf.mi_pval * ((n * n - n) / 2) * app_conf.ptest_size);

    // check if not too small
    if (null_size < 2) {
	std::cout << master << warning << "permutation test size too small"
		  << std::endl;

	// we want "second" largest value
	double r = app_conf.mi_pval * ((n * n - n) / 2);
	app_conf.ptest_size = static_cast<unsigned int>((2 / r) + 0.5);

	std::cout << master << "p-test changed, q = " << app_conf.ptest_size
		  << std::endl;
    }

    // get maximal size of local null dist
    data.r_lim
	= static_cast<unsigned int>(app_conf.mi_pval * app_conf.ptest_size * n * n);
} // check_test_size


bool run_bootstrap(const mpix::MPI_env& mpi_env, AppConfig& app_conf,
		   InputData& data, jaz::plain_array<double>& adj_matrix) {
    std::cout << master << ">>> bootstrap runs <<<" << std::endl;

    // set common seed
    long int seed = rand();
    MPI_Bcast(&seed, 1, MPI_LONG, mpi_env.root(), mpi_env.comm());
    srand(seed);

    jaz::plain_array<unsigned int> col_choice(data.g_ncol);
    jaz::plain_array<double> noise(data.g_ncol);
    jaz::numeric_table<double> exp_table(data.exp_table);

    if (app_conf.mem_report == true) MEM_REPORT;

    std::ostringstream ss;

    for (unsigned int i = 0; i < app_conf.boot_size; ++i) {
	std::cout << master << "--- run " << (i + 1) << "/"
		  << app_conf.boot_size << std::endl;

	// generate bootstrap data
	for (unsigned int j = 0; j < data.g_ncol; ++j) {
	    col_choice[j] = (rand() % data.g_ncol);
	    noise[j] = (9 - (rand() % 19)) * 1e-7;
	}

	for (unsigned int j = 0; j < data.num_rows; ++j) {
	    double* tab_in = exp_table.row_begin(j);
	    double* tab_out = data.exp_table.row_begin(j);
	    for (unsigned int k = 0; k < data.g_ncol; ++k) {
		tab_out[k] = tab_in[col_choice[k]];
		tab_out[k] = tab_in[col_choice[k]] + noise[k];
	    }

	}

	TIMER_START;

	// process new data
	rank_transform(mpi_env, app_conf, data);

	if (app_conf.mi_estim == AppConfig::BSPLINE) {
	    jaz::bs_mutual_info<double> I(data.g_ncol, app_conf.mi_b, app_conf.mi_k);
	    rw_build_adj_matrix(mpi_env, app_conf, data, adj_matrix, I);
	}
	else if (app_conf.mi_estim == AppConfig::GAUSSIAN) {
	    jaz::gk_mutual_info<double> I(data.g_ncol, jaz::h_silverman<double, 2>(data.g_ncol), jaz::h_silverman<double, 2>(data.g_ncol));
	    rw_build_adj_matrix(mpi_env, app_conf, data, adj_matrix, I);
	}

	rw_clean_adj_matrix(mpi_env, app_conf, data, adj_matrix);
	rw_post_adj_matrix(mpi_env, app_conf, data, adj_matrix);

	data.boot = i + 1;
	data.etime = static_cast<unsigned int>(TIMER_GET);

	ss.clear();
	ss.str("");
	ss << (i + 1);

	std::string name = app_conf.output + "." + ss.str();

	if (rw_write_adj(mpi_env, app_conf, data, adj_matrix, name) == false) {
	    return false;
	}

	// reset seed to be sync
	seed = rand();
	MPI_Bcast(&seed, 1, MPI_LONG, mpi_env.root(), mpi_env.comm());
	srand(seed);
    }

    return true;
} // run_bootstrap


bool run(const mpix::MPI_env& mpi_env, AppConfig& app_conf) {
    TIMER_START;

    app_conf.print_config(std::cout << master);

    if (app_conf.mem_report == true) MEM_REPORT;

    InputData data;
    if (app_conf.mi_pval == 1.0) data.has_null_dist = true;

    // read input data
    bool res = rw_read_exp(mpi_env, app_conf, data);
    if (res == false) return false;

    // sample if needed
    if (app_conf.sample_size > 0) sample_data(mpi_env, app_conf, data);

    if (app_conf.input_tf.empty() == false) {
    	res = read_tfs(mpi_env, app_conf, data);
    	if (res == false) return false;
	check_tfs(data);
    }

    if (app_conf.mem_report == true) MEM_REPORT;

    // rank transform data
    rank_transform(mpi_env, app_conf, data);

    if (app_conf.mem_report == true) MEM_REPORT;

    // check permutation test size
    check_test_size(app_conf, data);

    // compute adj matrix
    jaz::plain_array<double> adj_matrix;

    if (app_conf.input_mi.empty() == true) {

	unsigned int ncol = data.g_ncol;

	if (app_conf.mi_estim == AppConfig::BSPLINE) {
	    jaz::bs_mutual_info<double> I(ncol, app_conf.mi_b, app_conf.mi_k);
	    rw_build_adj_matrix(mpi_env, app_conf, data, adj_matrix, I);
	}
	else if (app_conf.mi_estim == AppConfig::GAUSSIAN) {
	    jaz::gk_mutual_info<double> I(ncol, jaz::h_silverman<double, 2>(ncol), jaz::h_silverman<double, 2>(ncol));
	    rw_build_adj_matrix(mpi_env, app_conf, data, adj_matrix, I);
	}

    } else {
	if (rw_read_adj(mpi_env, app_conf, data, adj_matrix) == false) {
	    return false;
	}
    }

    data.etime = static_cast<unsigned int>(TIMER_GET);
    if (app_conf.mem_report == true) MEM_REPORT;

    // remove edges between independent variables
    if (rw_clean_adj_matrix(mpi_env, app_conf, data, adj_matrix) == false) {
	return false;
    }

    // do post processing
    rw_post_adj_matrix(mpi_env, app_conf, data, adj_matrix);

    if (app_conf.mem_report == true) MEM_REPORT;

    // get time to create
    data.etime = static_cast<unsigned int>(TIMER_GET);

    // write final output
    res = rw_write_adj(mpi_env, app_conf, data, adj_matrix, app_conf.output);
    if (res == false) return false;

    if (app_conf.mem_report == true) MEM_REPORT;

    // write MCL labs
    if (app_conf.output_lab.empty() == false) {
	res = rw_write_lab(mpi_env, app_conf, data, adj_matrix, app_conf.output_lab);
	if (res == false) return false;
    }

    if (app_conf.boot_size > 0) {
	if (run_bootstrap(mpi_env, app_conf, data, adj_matrix) == false) {
	    return false;
	}
    }

    return true;
} // run


int main(int argc, char* argv[]) {
    const int MIN_CPU = 1;

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

    // HERE WE START
    std::cout << master << FULL_NAME << "\n"
	      << SHORT_NAME << " " << MAJOR_VERSION << "." << MINOR_VERSION
	      << " coded 2007-2010 by Jaroslaw Zola" << std::endl;

    AppConfig app_conf;

    if (app_conf.parse(argc, argv, std::cout << master, mpi_env.am_I_root()) == false) {
	return MPI_Finalize();
    }

    if (mpi_env.size() < MIN_CPU) {
	std::cout << master << error
		  << mpi_env.size() << " CPU available, "
		  << "more than " << MIN_CPU - 1 << " needed" << std::endl;
	return MPI_Finalize();
    }

    if (mpi_env.am_I_root() == true) {
	std::cout << "ARGV:";
	for (int i = 1; i < argc; ++i) std::cout << " " << argv[i];
	std::cout << ", CPU: " << mpi_env.size() << std::endl;
    }


    // MAIN PROCEDURE
    TIMER_START;

    srand(app_conf.rng_seed);

    // here we should call processor
    bool res = run(mpi_env, app_conf);

    if (res == true) {
	MPI_Barrier(MPI_COMM_WORLD);
	std::cout << master << "\nDone: " << timer(TIMER_GET) << "\n\n";
    } else {
	std::cout << every << "Stopped: " << timer(TIMER_GET) << "\n";
	exit(-1);
    }

    return MPI_Finalize();
} // main
