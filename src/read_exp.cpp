/***
 *  $Id: read_exp.cpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: read_exp.cpp
 *  Created: Dec 16, 2007
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

#include "read_exp.hpp"


bool create_index(const AppConfig& app_conf,
		  unsigned int& nrow, unsigned int& ncol, unsigned int& pvals,
		  std::vector<unsigned int>& index,
		  std::vector<std::string>& probes) {
    std::ifstream f(app_conf.input.c_str());

    if (!f) {
	std::cout << error << "can't open file " << app_conf.input << std::endl;
	return false;
    }

    // create custom buffer
    // can improve performance?
    // const unsigned int N = 8192 * 4096;
    // jaz::plain_array<char> fbuf(N);
    // f.rdbuf()->pubsetbuf(fbuf.begin(), N);

    nrow = 0;
    ncol = 0;

     // check header
    std::string s;
    std::getline(f, s);

    if (f.good() == false) {
	std::cout << error << "wrong header in input file" << std::endl;
	return false;
    }

    std::vector<std::string> tokens;
    std::back_insert_iterator<std::vector<std::string> > tokens_ii(tokens);

    jaz::split('\t', s, tokens_ii);
    ncol = tokens.size();

    unsigned int pos;

    std::string probe;
    std::istringstream ss;

    jaz::uc_compare uc_cmp;

    // check description
    do {
	pos = f.tellg();
	std::getline(f, s);

	ss.clear();
	ss.str(s);

	std::getline(ss, probe, '\t');

	if ((s.empty() == true) || (f.good() == false)) {
	    std::cout << error << "wrong header in input file" << std::endl;
	    return false;
	}
    } while (uc_cmp("Description", probe) == 0);

    std::set<std::string> pset;

    // this should be first line with data
    tokens.clear();
    jaz::split('\t', s, tokens_ii);
    pvals = (((tokens.size() - 2) >> 1) == (ncol - 2));

    // check main part
    do {
	++nrow;
	index.push_back(pos);

	ss.clear();
	ss.str(s);

	std::getline(ss, probe, '\t');

	if (pset.find(probe) != pset.end()) {
	    std::cout << error << "redundant entry " << probe
		      << " in data row " << nrow << std::endl;
	    return false;
	}
	else pset.insert(probe);

	probes.push_back(probe);

	pos = f.tellg();
	std::getline(f, s);
    }
    while ((f.good() == true) && (s.empty() == false));

    index.push_back(pos);
    ncol -= 2;

    f.close();

    return true;
} // create_index


bool create_table(const AppConfig& app_conf,
		  char* buf, unsigned int buf_sz, InputData& data, bool pvals) {
    std::istringstream ss;

    if (app_conf.has_pubsetbuf == true) ss.rdbuf()->pubsetbuf(buf, buf_sz);
    else ss.str(std::string(buf, buf_sz));

    unsigned int nrow = data.num_rows;
    unsigned int ncol = data.g_ncol;

    data.exp_table = jaz::numeric_table<double>(nrow, ncol);

    if (app_conf.mem_report == true) MEM_REPORT;

    std::string s;
    double* tab;

    double d;

    for (unsigned int i = 0; i < nrow; ++i) {
	std::getline(ss, s, '\t');
	std::getline(ss, s, '\t');

	if (ss.good() == false) {
	    std::cout << every << error << "incorrect value "
		      << s << " in data row "
		      << data.first_row + i + 1 << " ?" << std::endl;
	    return false;
	}

	tab = data.exp_table.row_begin(i);

	for (unsigned int j = 0; j < ncol; ++j) {
	    ss >> tab[j];

	    if (ss.good() == false) {
		std::cout << every << error << "incorrect value in line "
			  << data.first_row + i + 1 << " column "
			  << (pvals + 1) * j + 1 << std::endl;
		return false;
	    }

	    if (pvals == true) ss >> d;
	} // for ncol
    } // for nrow

    return true;
} // create_table


void send_probe_names(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
		      InputData& data) {
    unsigned int tot_len = data.g_nrow;

    if (mpi_env.am_I_root() == true) {
	for (unsigned int i = 0; i < data.g_nrow; ++i) {
	    tot_len += data.probes[i].size();
	}
    }

    MPI_Bcast(&tot_len, 1, MPI_UNSIGNED, mpi_env.root(), mpi_env.comm());

    char* buf = new char[tot_len];
    char* pos = buf;

    if (mpi_env.am_I_root() == true) {
	for (unsigned int i = 0; i < data.g_nrow; ++i) {
	    std::copy(data.probes[i].begin(), data.probes[i].end(), pos);
	    pos += data.probes[i].size();
	    *pos = '\t';
	    ++pos;
	}
    }

    MPI_Bcast(buf, tot_len, MPI_CHAR, mpi_env.root(), mpi_env.comm());

    if (mpi_env.am_I_root() == false) {
	data.probes.resize(data.g_nrow);

	std::istringstream ss;

	if (app_conf.has_pubsetbuf == true) ss.rdbuf()->pubsetbuf(buf, tot_len);
	else ss.str(std::string(buf, tot_len));

	for (unsigned int i = 0; i < data.g_nrow; ++i) {
	    std::getline(ss, data.probes[i], '\t');
	}
    }

    delete[] buf;
} // send_probe_names


bool rw_read_exp(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
		 InputData& data) {
    TIMER_START;
    std::cout << master << "* reading expression data..." << std::endl;

    unsigned int dim[3];
    std::vector<unsigned int> index;

    // get file layout
    if (mpi_env.am_I_root() == true) {
	char res = 1;

	if (create_index(app_conf, dim[0], dim[1], dim[2], index, data.probes) == false) {
	    res = 0;
	    MPI_Bcast(&res, 1, MPI_CHAR, mpi_env.root(), mpi_env.comm());
	    return false;
	}

	std::cout << "data size: " << dim[0] << " rows, " << dim[1]
		  << " columns";

	if (dim[2] != 0) std::cout << ", p-values included";
	std::cout << std::endl;

	if ((dim[0] < 2) || (dim[1] < 2)) {
	    std::cout << error << "data is too small" << std::endl;
	    res = 0;
	    MPI_Bcast(&res, 1, MPI_CHAR, mpi_env.root(), mpi_env.comm());
	    return false;
	}

	MPI_Bcast(&res, 1, MPI_CHAR, mpi_env.root(), mpi_env.comm());
    } else {
	char res = 1;
	MPI_Bcast(&res, 1, MPI_CHAR, mpi_env.root(), mpi_env.comm());
	if (res == 0) return false;
    } // if root


    // send dimension data
    MPI_Bcast(dim, 3, MPI_UNSIGNED, mpi_env.root(), mpi_env.comm());

    data.g_nrow = dim[0];
    data.g_ncol = dim[1];

    if (mpi_env.am_I_root() == false) index.resize(dim[0] + 1);

    MPI_Bcast(&index[0], dim[0] + 1, MPI_UNSIGNED, mpi_env.root(), mpi_env.comm());

    // send probes
    send_probe_names(mpi_env, app_conf, data);

    // read part of data
    unsigned int frow, lrow;
    data.num_rows = row2cpu(mpi_env.rank(), mpi_env.size(), dim[0], frow, lrow);

    data.first_row = frow;
    data.last_row = lrow;

    std::cout << master << "number of rows per cpu: ~" << data.num_rows
	      << std::endl;

    if (data.num_rows < 2) {
	std::cout << every << error << "too many CPUs for this task" << std::endl;
	return false;
    }

    unsigned int buf_sz = index[lrow] - index[frow];
    char* buf = new char[buf_sz];

    if (app_conf.mem_report == true) MEM_REPORT;

    MPI_File fh;
    MPI_Status mpi_stat;

    MPI_File_open(mpi_env.comm(), const_cast<char*>(app_conf.input.c_str()),
		  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    MPI_File_seek(fh, index[frow], MPI_SEEK_SET);
    MPI_File_read(fh, buf, buf_sz, MPI_CHAR, &mpi_stat);

    MPI_File_close(&fh);

    // parse buffered data
    if (create_table(app_conf, buf, buf_sz, data, dim[2]) == false) return false;

    delete[] buf;

    std::cout << master << "reading done: " << timer(TIMER_GET) << std::endl;

    return true;
} // rw_read_exp
