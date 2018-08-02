/***
 *  $Id: adj_build.hpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: adj_build.hpp
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

#ifndef ADJ_BUILD_HPP
#define ADJ_BUILD_HPP

#include <mpix/MPI_env.hpp>
#include <jaz/mutual_info.hpp>
#include <jaz/stat_util.hpp>
#include <jaz/utility_add.hpp>
#include "AppConfig.hpp"
#include "InputData.hpp"
#include "mpi_tags.hpp"
#include "write_adj.hpp"
#include "sort_reduce.hpp"
#include "utility.hpp"
#include <string.h>


void exchange_block(const mpix::MPI_env&, const AppConfig&, const InputData&,
		    jaz::plain_array<double>&);

std::pair<double, double> find_value_range(const mpix::MPI_env&,
					   const jaz::plain_array<double>&);

unsigned int get_tot_edges(const mpix::MPI_env&, const InputData&,
			   const jaz::plain_array<double>&);

bool rw_clean_adj_matrix(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
			 InputData& data, jaz::plain_array<double>& adj_matrix);

void rw_post_adj_matrix(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
			InputData& data, jaz::plain_array<double>& adj_matrix);


// structure to describe block of input data
// received from remote processor
template <typename T> struct DataBlock {
    unsigned int frow;
    unsigned int lrow;
    unsigned int num_rows;

    jaz::plain_array<T> data;
}; // struct DataBlock


template <typename MIEstim> class p_test {
public:
    p_test(MIEstim* mi, unsigned int n, unsigned int k)
	: I_(mi), n_(n), k_(k), yi_(n) { }

    void operator()(const unsigned int* yi, InputData& data) {
	memcpy(yi_.begin(), yi, n_ * sizeof(unsigned int));

	for (unsigned int i = 0; i < k_; ++i) {
	    std::random_shuffle(yi_.begin(), yi_.end());
	    double val = std::max(0.0, (*I_)(yi_.begin()));
	    data.null_dist.push_back(val);
	    data.nd_sum += val;
	    data.nd_sum_sqr += (val * val);
	    ++data.nd_S;
	}

	std::sort(data.null_dist.begin(), data.null_dist.end(),
		  std::greater<double>());

	if (data.r_lim < data.null_dist.size()) {
	    data.null_dist.resize(data.r_lim);
	}
    } // operator()

private:
    MIEstim* I_;

    unsigned int n_;
    unsigned int k_;

    jaz::plain_array<unsigned int> yi_;

}; // class p_test


template <typename MIEstim>
void compute_local_block(const AppConfig& app_conf,
			 InputData& data, jaz::plain_array<double>& adj_matrix,
			 MIEstim& I, p_test<MIEstim>& pt) {
    unsigned int n = data.num_rows;

    unsigned int g_row = data.first_row;
    unsigned int g_col = 0;

    unsigned int ncol = data.g_ncol;

    for (unsigned int i = 0; i < n - 1; ++i, ++g_row) {
	g_col = g_row + 1;

	adj_matrix[i * data.g_nrow + g_col]
	    = adj_matrix[(i + 1) * data.g_nrow + g_row]
	    = I(&data.rank_array[i * ncol], &data.rank_array[(i + 1) * ncol]);

	// ptest
	if (data.has_null_dist == false) {
	    pt(&data.rank_array[(i + 1) * ncol], data);
	}

	++g_col;

	for (unsigned int j = i + 2; j < n; ++j, ++g_col) {
	    adj_matrix[i * data.g_nrow + g_col]
		= adj_matrix[j * data.g_nrow + g_row]
		= I(&data.rank_array[j * ncol]);

	    // ptest
	    if (data.has_null_dist == false) {
		pt(&data.rank_array[j * ncol], data);
	    }
	}
    }
} // compute_local_block


template <typename MIEstim>
void compute_block(const AppConfig& app_conf,
		   InputData& data, const DataBlock<unsigned int>& recv_data,
		   jaz::plain_array<double>& adj_matrix,
		   MIEstim& I, p_test<MIEstim>& pt) {
    unsigned int n = data.num_rows;
    unsigned int m = recv_data.num_rows;

    unsigned int ncol = data.g_ncol;

    for (unsigned int i = 0; i < n; ++i) {
	unsigned int g_col = recv_data.frow;

	adj_matrix[i * data.g_nrow + g_col] =
	    I(&data.rank_array[i * ncol], &recv_data.data[0]);

	// ptest
	if (data.has_null_dist == false) pt(&recv_data.data[0], data);

	++g_col;

	for (unsigned int j = 1; j < m; ++j, ++g_col) {
	    adj_matrix[i * data.g_nrow + g_col] = I(&recv_data.data[j * ncol]);

	    // ptest
	    if (data.has_null_dist == false) pt(&recv_data.data[j * ncol], data);

	}
    }
} // compute_block


template <typename MIEstim>
void compute_half_block(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
			InputData& data, const DataBlock<unsigned int>& recv_data,
			jaz::plain_array<double>& adj_matrix,
			MIEstim& I, p_test<MIEstim>& pt) {
    unsigned int n = data.num_rows;
    unsigned int m = recv_data.num_rows;

    unsigned int ncol = data.g_ncol;

    unsigned int n_min = 0;
    unsigned int n_max = n;
    unsigned int m_min = 0;
    unsigned int m_max = m;

    if (mpi_env.rank() < (mpi_env.size() >> 1)) {
	// to process half-block above diagonal
	n_max = static_cast<unsigned int>(((1.0 * n) / 2) + 0.5);
    } else {
	// process half-block below diagonal
	m_min = m - (m / 2);
    }

    // here goes processing
    for (unsigned int i = n_min; i < n_max; ++i) {
	unsigned int g_col = recv_data.frow + m_min;

	adj_matrix[i * data.g_nrow + g_col] =
	    I(&data.rank_array[i * ncol], &recv_data.data[m_min * ncol]);

	// ptest
	if (data.has_null_dist == false) pt(&recv_data.data[m_min * ncol], data);

	++g_col;

	for (unsigned int j = m_min + 1; j < m_max; ++j, ++g_col) {
	    adj_matrix[i * data.g_nrow + g_col] = I(&recv_data.data[j * ncol]);

	    // ptest
	    if (data.has_null_dist == false) pt(&recv_data.data[j * ncol], data);

	}
    }

} // compute_half_block


template <typename MIEstim>
void rw_build_adj_matrix(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
			 InputData& data, jaz::plain_array<double>& adj_matrix,
			 MIEstim& I) {
    std::cout << master << "* building adjacency matrix..." << std::endl;
    TIMER_START;

    unsigned int nrow = data.g_nrow;
    unsigned int ncol = data.g_ncol;

    // we use min double to mark empty cells
    // we use this value frequently later
    adj_matrix.resize(data.num_rows * nrow);
    for (unsigned int i = 0; i < adj_matrix.size(); i++) {
	adj_matrix[i] = std::numeric_limits<double>::min();
    }


    // *** and p_tester ***
    p_test<MIEstim> pt(&I, ncol, app_conf.ptest_size);


    if (app_conf.mem_report == true) MEM_REPORT;

    MPI_Status mpi_stat;

    int rank = mpi_env.rank();
    int size = mpi_env.size();

    unsigned int to = (size + rank - 1) % size;
    unsigned int from = (rank + 1) % size;

    DataBlock<unsigned int> recv_data;
    DataBlock<unsigned int> send_data;


    // COMPUTE BLOCKS OF ADJ_MATRIX
    compute_local_block(app_conf, data, adj_matrix, I, pt);

    if (mpi_env.size() > 1) {
	// make first shift
	double t0 = MPI_Wtime();

	recv_data.num_rows = row2cpu(from, size, nrow,
				     recv_data.frow, recv_data.lrow);

	unsigned int recv_size = recv_data.num_rows * ncol;

	recv_data.data.resize(recv_size);

	double bsT0 = MPI_Wtime();

	MPI_Sendrecv(data.rank_array.begin(), data.num_rows * ncol, MPI_UNSIGNED,
		     to, ADJ_BUILD_DATA_BLOCK,
		     recv_data.data.begin(), recv_size, MPI_UNSIGNED,
		     from, ADJ_BUILD_DATA_BLOCK, mpi_env.comm(), &mpi_stat);

	double bsT1 = MPI_Wtime();

	if (mpi_env.size() != 2) {
	    compute_block(app_conf, data, recv_data, adj_matrix, I, pt);
	} else compute_half_block(mpi_env, app_conf, data, recv_data, adj_matrix, I, pt);

	// we can avoid extra work for odd cases
	// in even cases we will have one iteration extra
	// thus each cpu will process only half of the block
	unsigned int lim = (mpi_env.size() >> 1);
	bool has_hblck = !(mpi_env.size() % 2);

	double T = (lim - 1) * (MPI_Wtime() - t0);

	if (app_conf.mem_report == true) MEM_REPORT;

	std::cout << master << "block shift time: "
		  << timer(bsT1 - bsT0) << std::endl;

	std::cout << master << "estimated remaining time: "
		  << timer(std::max(T, 0.0)) << std::endl;

	// proceed with shifting
	DataBlock<unsigned int>* sblock = &recv_data;
	DataBlock<unsigned int>* rblock = &send_data;

	// be verbose if needed
	unsigned int v_step
	    = static_cast<unsigned int>(ceil(static_cast<double>(lim) / 10));
	if (app_conf.verbose == true) std::cout << master << "." << std::flush;

	for (unsigned int i = 1; i < lim; ++i) {
	    unsigned int send_size = recv_size;

	    // get recv_size
	    unsigned int data_rank = (rank + 1 + i) % size;
	    (*rblock).num_rows = row2cpu(data_rank, size, nrow,
					 (*rblock).frow, (*rblock).lrow);

	    recv_size = (*rblock).num_rows * ncol;
	    (*rblock).data.resize(recv_size);

	    // shift
	    MPI_Sendrecv((*sblock).data.begin(), send_size, MPI_UNSIGNED,
			 to, ADJ_BUILD_DATA_BLOCK,
			 (*rblock).data.begin(), recv_size, MPI_UNSIGNED,
			 from, ADJ_BUILD_DATA_BLOCK, mpi_env.comm(), &mpi_stat);

	    if (has_hblck && (i == lim - 1)) {
		compute_half_block(mpi_env, app_conf, data, *rblock, adj_matrix, I, pt);
	    } else compute_block(app_conf, data, *rblock, adj_matrix, I, pt);

	    std::swap(sblock, rblock);

	    if (app_conf.verbose == true) {
		if ((i % v_step) == 0) std::cout << master << "." << std::flush;
	    }
	} // for i

	// EXCHANGE BLOCKS TO GET COMPLETE ADJ_MATRIX ROW WISE DISTRIBUTED
	exchange_block(mpi_env, app_conf, data, adj_matrix);

    } // if (mpi_env.size() > 1)


    // DEBUG JUNK
    // std::cout << std::setprecision(3) << std::setfill(' ') << std::left;
    // for (unsigned int i = 0; i < mpi_env.size(); ++i) {
    // 	if (mpi_env.rank() == i) {
    // 	    for (unsigned int x = 0; x < data.num_rows; ++x) {
    // 		for (unsigned int y = 0; y < data.g_nrow; ++y) {
    // 		    std::cout << std::setw(9) << adj_matrix[x * data.g_nrow + y] << " ";
    // 		}
    // 		std::cout << std::endl;
    // 	    }
    // 	    std::cout << "----------" << std::endl;
    // 	}
    // 	MPI_Barrier(mpi_env.comm());
    // } // for i

    // find min and max value
    std::pair<double, double> g_val = find_value_range(mpi_env, adj_matrix);

    std::cout << master << "values range: [" << g_val.first << ","
	      << g_val.second << "]" << std::endl;

    data.has_null_dist = true;

    std::cout << master << "building done: " << timer(TIMER_GET) << std::endl;
} // rw_build_adj_matrix


#endif // ADJ_BUILD_HPP
