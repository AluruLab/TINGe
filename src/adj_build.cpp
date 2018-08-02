/***
 *  $Id: adj_build.cpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: adj_build.cpp
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

#include "adj_build.hpp"


// to use template in MPI_Op
extern "C" void dg_sort_reduce(void* in, void* inout, int* len, MPI_Datatype* d) {
    sort_reduce<double, std::greater>(in, inout, len, d);
}


void exchange_block(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
		    const InputData& data, jaz::plain_array<double>& adj_matrix) {
    TIMER_START;

    int rank = mpi_env.rank();
    int size = mpi_env.size();

    MPI_Status mpi_stat;

    unsigned int nrow = data.num_rows;

    jaz::plain_array<double> send_buf;
    jaz::plain_array<double> recv_buf;

    unsigned int lim = (size >> 1) + (size % 2);

    // verbosity step
    unsigned int v_step
	= static_cast<unsigned int>(ceil(static_cast<double>(lim) / 10));
    if (app_conf.verbose == true) std::cout << master << "," << std::flush;

    for (unsigned int p = 0; p < lim; ++p) {
	unsigned int to = (rank + 1 + p) % size;
	unsigned int from = (size + rank - p - 1) % size;

	unsigned int to_frow, to_lrow;
	unsigned int to_nrow = row2cpu(to, size, data.g_nrow, to_frow, to_lrow);

	unsigned int send_size = nrow * to_nrow;
	send_buf.resize(send_size);

	for (unsigned int i = 0; i < nrow; ++i) {
	    for (unsigned int j = 0; j < to_nrow; ++j) {
		send_buf[i * to_nrow + j] = adj_matrix[i * data.g_nrow + to_frow + j];
	    }
	}

	unsigned int from_frow, from_lrow;
	unsigned int from_nrow = row2cpu(from, size, data.g_nrow,
				     from_frow, from_lrow);

	unsigned int recv_size = nrow * from_nrow;
	recv_buf.resize(recv_size);

	MPI_Sendrecv(send_buf.begin(), send_size, MPI_DOUBLE,
		     to, ADJ_BUILD_DATA_BLOCK,
		     recv_buf.begin(), recv_size, MPI_DOUBLE,
		     from, ADJ_BUILD_DATA_BLOCK,
		     mpi_env.comm(), &mpi_stat);

	// transpose and store
	for (unsigned int i = 0; i < nrow; ++i) {
	    for (unsigned int j = 0; j < from_nrow; ++j) {
		if (recv_buf[j * nrow + i] > std::numeric_limits<double>::min()) {
		    adj_matrix[i * data.g_nrow + from_frow + j] = recv_buf[j * nrow + i];
		}
	    }
	}

	// verbose output
	if (app_conf.verbose == true) {
	    if ((p % v_step) == 0) std::cout << master << "." << std::flush;
	}
    } // for p
    if (app_conf.verbose == true) std::cout << master << "\n";

    std::cout << master << "block exchange done: " << timer(TIMER_GET) << std::endl;
} // exchange_block


std::pair<double, double>
find_value_range(const mpix::MPI_env& mpi_env,
		 const jaz::plain_array<double>& adj_matrix) {

    std::pair<const double*, const double*> min_max =
	min_max_val(adj_matrix.begin(), adj_matrix.end(),
		    std::numeric_limits<double>::min());

    std::pair<double, double> g_val;

    MPI_Reduce(const_cast<double*>(min_max.first), &g_val.first, 1,
	       MPI_DOUBLE, MPI_MIN, mpi_env.root(), mpi_env.comm());

    MPI_Reduce(const_cast<double*>(min_max.second), &g_val.second, 1,
	       MPI_DOUBLE, MPI_MAX, mpi_env.root(), mpi_env.comm());

    return g_val;
} // find_value_range


unsigned int get_tot_edges(const mpix::MPI_env& mpi_env, const InputData& data,
			   const jaz::plain_array<double>& adj_matrix) {
    unsigned int tot_edg = 0;
    unsigned int num_edg = 0;

    for (unsigned int i = 0; i < data.num_rows; ++i) {
	for (unsigned int j = data.first_row + i + 1; j < data.g_nrow; ++j) {
	    if (adj_matrix[i * data.g_nrow + j] > 0.0) ++num_edg;
	}
    }

    MPI_Reduce(&num_edg, &tot_edg, 1, MPI_UNSIGNED, MPI_SUM,
	       mpi_env.root(), mpi_env.comm());

    return tot_edg;
} // get_tot_edges


std::pair<double, double>
find_mi_stats(const mpix::MPI_env& mpi_env, const InputData& data,
	      const jaz::plain_array<double>& adj_matrix, unsigned int tot_edg) {
    double x[2] = { 0.0, 0.0 };

    for (unsigned int i = 0; i < data.num_rows; ++i) {
	for (unsigned int j = data.first_row + i + 1; j < data.g_nrow; ++j) {
	    double v = adj_matrix[i * data.g_nrow + j];
	    if (v > 0.0) {
		x[0] += v;
		x[1] += v * v;
	    }
	}
    }

    double X[2] = { 0.0, 0.0 };
    MPI_Allreduce(&x, &X, 2, MPI_DOUBLE, MPI_SUM, mpi_env.comm());

    double mean = X[0] / tot_edg;
    double sd = sqrt((X[1] - 2 * mean * X[0] + tot_edg * mean * mean) / (tot_edg - 1));

    return std::pair<double, double>(mean, sd);
} // find_mi_stats



double find_pval_threshold(const mpix::MPI_env& mpi_env,
			   const AppConfig& app_conf,
			   InputData& data,
			   jaz::plain_array<double>& adj_matrix) {
    TIMER_START;

    unsigned int S = 0;
    unsigned int n = data.nd_S;

    if (n == 0) return 0.0;

    MPI_Allreduce(&n, &S, 1, MPI_UNSIGNED, MPI_SUM, mpi_env.comm());

    // compute mean and standard deviation of null distribution
    double nd_sums[2] = { data.nd_sum, data.nd_sum_sqr };
    double nd_sums_all[2] = { 0.0, 0.0 };

    MPI_Allreduce(nd_sums, nd_sums_all, 2, MPI_DOUBLE, MPI_SUM, mpi_env.comm());

    double mean = data.nd_mean = (nd_sums_all[0] / S);
    data.nd_sd = sqrt((nd_sums_all[1] - 2 * mean * nd_sums_all[0] + S * mean * mean)
		      / (S - 1));

    std::cout << master << "sample size: " << S << std::endl;
    std::cout << master << "random MI, mean: " << data.nd_mean << std::endl;
    std::cout << master << "random MI, stdev: " << data.nd_sd << std::endl;

    unsigned int r = static_cast<unsigned int>(app_conf.mi_pval * S);

    jaz::plain_array<double> buf;

    buf.push_back(data.null_dist.begin(), std::min<std::size_t>(r, data.null_dist.size()));
    buf.resize(r, 0.0);

    MPI_Op op;
    MPI_Op_create(dg_sort_reduce, 1, &op);

    jaz::plain_array<double> recv_buf(r);

    MPI_Reduce(buf.begin(), recv_buf.begin(), r, MPI_DOUBLE, op,
	       mpi_env.root(), mpi_env.comm());

    MPI_Op_free(&op);

    double thr = *(recv_buf.end() - 1);
    data.mi_thr_used = thr;

    MPI_Bcast(&thr, 1, MPI_DOUBLE, mpi_env.root(), mpi_env.comm());

    std::cout << master << "random MI, max: " << *(recv_buf.begin()) << std::endl;
    std::cout << master << "threshold value: " << thr << ", r = " << r << std::endl;

    std::cout << master << "threshold done: " << timer(TIMER_GET) << std::endl;

    return thr;
} // find_pval_threshold


void analyze_dpi_rows(const AppConfig& app_conf, const InputData& data,
		      unsigned int x, unsigned int z,
		      const double* row_x, const double* row_z,
		      std::vector<bool>& dmap) {

    unsigned int n = data.g_nrow;
    unsigned int lx = x - data.first_row;

    double i_xz = (row_x[z]) / (1.0 - app_conf.dpi_tol);

    for (unsigned int y = 0; y < n; ++y) {
	double i_xy = row_x[y];
	if (i_xy < i_xz) continue;

	double i_yz = row_z[y];
	if (i_yz < i_xz) continue;

	// check TFs constraints
	if (data.has_TF == true) {
	    if (data.isTF[y] == false) {
		if (((data.isTF[x] && data.isTF[z]) == false) &&
		    ((data.isTF[x] || data.isTF[z]) == true)) continue;
	    }
	}

	dmap[lx * n + z] = true;
	break;
    }
} // analyze_dpi_rows


void analyze_dpi_block(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
		       const InputData& data,
		       const jaz::plain_array<double>& adj_matrix,
		       std::vector<bool>& dmap) {
    int rank = mpi_env.rank();
    int size = mpi_env.size();

    MPI_Status mpi_stat;

    unsigned int nrow = data.g_nrow;

    int to = rank - 1;
    int from = rank + 1;

    DataBlock<double> recv_data;

    recv_data.num_rows = row2cpu(from, size, nrow, recv_data.frow, recv_data.lrow);

    unsigned int recv_size = recv_data.num_rows * nrow;

    // send first block
    if (rank == 0) {

	recv_data.data.resize(recv_size);
	if (app_conf.mem_report == true) MEM_REPORT;

	MPI_Recv(recv_data.data.begin(), recv_size, MPI_DOUBLE,
		 from, ADJ_BUILD_DATA_BLOCK, mpi_env.comm(), &mpi_stat);

    } else if (rank == (size - 1)) {

	MPI_Send(const_cast<double*>(adj_matrix.begin()), adj_matrix.size(),
		 MPI_DOUBLE, to, ADJ_BUILD_DATA_BLOCK, mpi_env.comm());

    } else {

	recv_data.data.resize(recv_size);

	MPI_Sendrecv(const_cast<double*>(adj_matrix.begin()), adj_matrix.size(),
		     MPI_DOUBLE, to, ADJ_BUILD_DATA_BLOCK,
		     recv_data.data.begin(), recv_size,
		     MPI_DOUBLE, from, ADJ_BUILD_DATA_BLOCK,
		     mpi_env.comm(), &mpi_stat);

    } // if


    if (rank < (size - 1)) {
	for (unsigned int x = 0; x < data.num_rows; ++x) {
	    const double* row_x = adj_matrix.begin() + x * nrow;

	    for (unsigned int z = 0; z < recv_data.num_rows; ++z) {
		const double* row_z = recv_data.data.begin() + z * nrow;

		if (row_x[recv_data.frow + z] > 0.0) {
		    analyze_dpi_rows(app_conf, data,
				     x + data.first_row, recv_data.frow + z,
				     row_x, row_z, dmap);
		}
	    }
	}
    }


    DataBlock<double> send_data;

    DataBlock<double>* rblock = &send_data;
    DataBlock<double>* sblock = &recv_data;

    // proceed with remaining data
    for (int i = 1; i < size - 1; ++i) {
	if (rank == 0) rblock = &recv_data;

	(*rblock).num_rows = row2cpu(from + i, size, nrow,
				     (*rblock).frow, (*rblock).lrow);

	recv_size = (*rblock).num_rows * nrow;

	if (rank == 0) {
	    recv_data.data.resize(recv_size);

	    MPI_Recv(recv_data.data.begin(), recv_size, MPI_DOUBLE,
		     from, ADJ_BUILD_DATA_BLOCK, mpi_env.comm(), &mpi_stat);

	} else if (rank == (size - i - 1)) {

	    MPI_Send((*sblock).data.begin(), (*sblock).data.size(), MPI_DOUBLE,
		     to, ADJ_BUILD_DATA_BLOCK, mpi_env.comm());

	} else if (rank < (size - i - 1)) {
	    (*rblock).data.resize(recv_size);

	    MPI_Sendrecv((*sblock).data.begin(), (*sblock).data.size(), MPI_DOUBLE,
			 to, ADJ_BUILD_DATA_BLOCK,
			 (*rblock).data.begin(), recv_size, MPI_DOUBLE,
			 from, ADJ_BUILD_DATA_BLOCK,
			 mpi_env.comm(), &mpi_stat);
	}

	// process new data
	if (rank < (size - i - 1)) {
	    for (unsigned int x = 0; x < data.num_rows; ++x) {
		const double* row_x = adj_matrix.begin() + x * nrow;

		for (unsigned int z = 0; z < (*rblock).num_rows; ++z) {
		    const double* row_z = (*rblock).data.begin() + z * nrow;

		    if (row_x[(*rblock).frow + z] > 0.0) {
			analyze_dpi_rows(app_conf, data,
					 x + data.first_row, (*rblock).frow + z,
					 row_x, row_z, dmap);
		    }
		}
	    }
	}

	std::swap(rblock, sblock);
    } // for mpi_env.size()
} // analyze_dpi_block


void apply_dpi(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
	       const InputData& data, jaz::plain_array<double>& adj_matrix) {
    TIMER_START;

    std::cout << master << "* applying DPI..." << std::endl;

    unsigned int nrow = data.num_rows;
    unsigned int ncol = data.g_nrow;

    std::vector<bool> dmap(nrow * ncol, false);

    // process local block
    for (unsigned int x = 0; x < nrow - 1; ++x) {
	const double* row_x = adj_matrix.begin() + x * ncol;

	for (unsigned int z = x + 1; z < nrow; ++z) {
	    const double* row_z = adj_matrix.begin() + z * ncol;

	    if (row_x[z + data.first_row] > 0.0) {
		analyze_dpi_rows(app_conf, data,
				 x + data.first_row, z + data.first_row,
				 row_x, row_z, dmap);
	    }
	}
    }

    // block process
    if (mpi_env.size() > 1) {
	analyze_dpi_block(mpi_env, app_conf, data, adj_matrix, dmap);
    }

    // make cleaning
    for (unsigned int x = 0; x < nrow; ++x) {
	for (unsigned int z = data.first_row + x + 1; z < ncol; ++z) {
	    if (dmap[x * ncol + z] == true) adj_matrix[x * ncol + z] = 0.0;
	}
    }

    std::cout << master << "DPI done: " << timer(TIMER_GET) << std::endl;
} // apply_dpi


bool rw_clean_adj_matrix(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
			 InputData& data, jaz::plain_array<double>& adj_matrix) {
    std::cout << master << "* cleaning adjacency matrix..." << std::endl;

    double mi_thr = app_conf.mi_thr;

    if (app_conf.mi_pval < 1.0) {
	mi_thr = find_pval_threshold(mpi_env, app_conf, data, adj_matrix);
    }

    // remove edges less than threshold
    for (unsigned int i = 0; i < data.num_rows; ++i) {
	for (unsigned int j = 0; j < data.g_nrow; ++j) {
	    if (adj_matrix[i * data.g_nrow + j] < mi_thr) {
		adj_matrix[i * data.g_nrow + j] = 0.0;
	    }
	}
    }

    // get network properties
    unsigned int max_edg = ((data.g_nrow * data.g_nrow) - data.g_nrow);
    unsigned int tot_edg = get_tot_edges(mpi_env, data, adj_matrix);

    // get density
    if (mpi_env.am_I_root()) {
	double d = (2.0 * tot_edg) / max_edg;
	double fdr = (app_conf.mi_pval * max_edg) / (2.0 * tot_edg);
	std::cout << "network density: " << d << std::endl;
	std::cout << "estimated FDR: " << fdr << std::endl;
	data.fdr = fdr;
    }

    // write relevance network
    if (app_conf.output_mi.empty() == false) {
	jaz::tie(data.mi_mean, data.mi_sd) = find_mi_stats(mpi_env, data, adj_matrix, tot_edg);
	bool res = rw_write_adj(mpi_env, app_conf, data, adj_matrix, app_conf.output_mi);
	if (res == false) return false;
    }

    // apply DPI
    if (app_conf.dpi_tol < 1.0) {
	apply_dpi(mpi_env, app_conf, data, adj_matrix);
	tot_edg = get_tot_edges(mpi_env, data, adj_matrix);

	// get density
	if (mpi_env.am_I_root()) {
	    double d = (2.0 * tot_edg) / max_edg;
	    std::cout << "new network density: " << d << std::endl;
	}

	// to compute final FDR when writing to file
	data.fdr = -1.0;
    } // DPI

    // get stats of MI
    jaz::tie(data.mi_mean, data.mi_sd) = find_mi_stats(mpi_env, data, adj_matrix, tot_edg);
    std::cout << master << "final MI, mean: " << data.mi_mean << std::endl;
    std::cout << master << "final MI, stdev: " << data.mi_sd << std::endl;

    std::cout << master << "cleaning done" << std::endl;
    return true;
} // rw_clean_adj_matrix


void rw_post_adj_matrix(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
			InputData& data, jaz::plain_array<double>& adj_matrix) {
    if (app_conf.mi_conv == AppConfig::CORRELATION) {
	std::cout << master << "* converting to correlation coefficients" << std::endl;
	for (unsigned int i = 0; i < adj_matrix.size(); ++i) {
	    adj_matrix[i] = sqrt(1.0 - exp(-2 * adj_matrix[i]));
	}
    }
} // rw_post_adj_matrix
