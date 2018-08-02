/***
 *  $Id: read_adj.cpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: read_adj.cpp
 *  Created: Apr 28, 2008
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

#include "read_adj.hpp"


typedef std::map<std::string, unsigned int> smap_type;
// typedef std::tr1::unordered_map<std::string, unsigned int> smap_type;


struct adj_block {
    double val;
    int id;
}; // struct adj_block


template <typename Map>
bool process_block(std::istream& is, const Map& probe2id,
		   jaz::plain_array<adj_block>& buf) {
    char c;
    std::string s;

    smap_type::const_iterator probe_iter;

    unsigned int id1, id2;
    double val;

    std::getline(is, s, '\t');
    probe_iter = probe2id.find(s);

    if (probe_iter == probe2id.end()) {
	std::cout << error << "unknown entry " << s << std::endl;
	return false;
    } else id1 = probe_iter->second;

    buf[0].id = id1;
    buf[0].val = -1.0;

    unsigned int pos = 1;

    while (is.good() == true) {
	std::getline(is, s, '\t');
	is >> val;

	if (!is) {
	    std::cout << error << "wrong value for edge " << s << std::endl;
	    return false;
	}

	is.get(c);

	probe_iter = probe2id.find(s);

	if (probe_iter == probe2id.end()) {
	    std::cout << every << error << "uknown entry " << s << std::endl;
	    return false;
	} else id2  = probe_iter->second;

	if (buf.size() < pos + 1) buf.resize(buf.size() + 1);

	buf[pos].id = id2;
	buf[pos].val = val;

	++pos;
    } // while is

    buf.resize(pos);

    return true;
} // process_block


void store_buffer(const InputData& data, const jaz::plain_array<adj_block>& buf,
		  jaz::plain_array<double>& adj_matrix) {
    unsigned int nrow = data.g_nrow;

    unsigned int id1 = buf[0].id;
    unsigned int id2;

    for (unsigned int i = 1; i < buf.size(); ++i) {
	id2 = buf[i].id;

	if (buf[i].val < 0.0) {
	    id1 = id2;
	    continue;
	}

	if ((data.first_row <= id1) && (id1 < data.last_row)) {
	    adj_matrix[(id1 - data.first_row) * nrow + id2] = buf[i].val;
	}

	if ((data.first_row <= id2) && (id2 < data.last_row)) {
	    adj_matrix[(id2 - data.first_row) * nrow + id1] = buf[i].val;
	}
    }

} // store_buffer


bool rw_read_adj(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
		 const InputData& data, jaz::plain_array<double>& adj_matrix) {
    TIMER_START;
    std::cout << master << "* reading adjacency matrix..." << std::endl;

    int rank = mpi_env.rank();
    jaz::plain_array<unsigned int> index;

    // generate index
    if (mpi_env.am_I_root() == true) {
	std::ifstream f(app_conf.input_mi.c_str());

	if (!f) {
	    std::cout << error << "can't open file "
		      << app_conf.input_mi << std::endl;
	    return false;
	}

	std::string s;
	unsigned int pos;

	do {
	    pos = f.tellg();
	    std::getline(f, s);
	}
	while ((f.eof() == false) && (s.empty() == false) && (s[0] == '>'));

	if (f.eof() == true) {
	    std::cout << error << "wrong file format" << std::endl;
	    return false;
	}

	do {
	    index.push_back(pos);
	    pos = f.tellg();
	    std::getline(f, s);
	}
	while (f.eof() == false);

	index.push_back(pos);

	f.close();
    } // root

    // bcast index
    unsigned int index_sz = index.size();
    MPI_Bcast(&index_sz, 1, MPI_UNSIGNED, mpi_env.root(), mpi_env.comm());

    index.resize(index_sz);
    MPI_Bcast(index.begin(), index_sz, MPI_UNSIGNED, mpi_env.root(), mpi_env.comm());

    // get block to read
    unsigned int frow, lrow;
    unsigned int nrow = row2cpu(rank, mpi_env.size(), index_sz - 1, frow, lrow);

    // read block
    MPI_File fh;
    MPI_Status mpi_stat;

    unsigned int buf_sz = 0;
    jaz::plain_array<char> buf;

    if (nrow > 0) buf_sz = index[lrow] - index[frow];
    buf.resize(buf_sz, 0);

    if (app_conf.mem_report == true) MEM_REPORT;

    MPI_File_open(mpi_env.comm(), const_cast<char*>(app_conf.input_mi.c_str()),
		  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    MPI_File_seek(fh, index[frow], MPI_SEEK_SET);
    MPI_File_read(fh, buf.begin(), buf_sz, MPI_CHAR, &mpi_stat);

    MPI_File_close(&fh);

    // process data
    adj_matrix.resize(data.num_rows * data.g_nrow, 0.0);
    jaz::plain_array<adj_block> dbuf(data.g_nrow);

    std::map<std::string, unsigned int> probe2id;
    for (unsigned int i = 0; i < data.g_nrow; ++i) probe2id[data.probes[i]] = i;

    if (app_conf.mem_report == true) MEM_REPORT;

    std::string s;
    std::istringstream is;
    std::istringstream is_b;

    if (app_conf.has_pubsetbuf == true) is.rdbuf()->pubsetbuf(buf.begin(), buf_sz);
    else is.str(std::string(buf.begin(), buf_sz));

    // get number of communication rounds
    unsigned int comm_rnd = row2cpu(0, mpi_env.size(), index_sz - 1, frow, lrow);

    jaz::plain_array<adj_block> rbuf;
    jaz::plain_array<int> dbuf_sz(mpi_env.size());
    jaz::plain_array<int> dbuf_disp(mpi_env.size());

    for (unsigned int i = 0; i < comm_rnd; ++i) {
	std::getline(is, s);

	dbuf_sz[rank] = 0;

	if (is.good() == true) {
	    is_b.clear();
	    is_b.str(s);

	    // process line
	    process_block(is_b, probe2id, dbuf);

	    dbuf_sz[rank] = dbuf.size();
	}

	MPI_Allgather(dbuf_sz.begin() + rank, 1, MPI_INT,
		      dbuf_sz.begin(), 1, MPI_INT, mpi_env.comm());

	unsigned int len = 0;
	dbuf_disp[0] = 0;

	for (int j = 1; j < mpi_env.size(); ++j) {
	    len += dbuf_sz[j - 1];
	    dbuf_disp[j] = len;
	}

	len += dbuf_sz[mpi_env.size() - 1];
	rbuf.resize(len);

	if (dbuf_sz[rank] > 0) {
	    MPI_Allgatherv(dbuf.begin(), dbuf.size(), MPI_DOUBLE_INT,
			   rbuf.begin(), dbuf_sz.begin(), dbuf_disp.begin(),
			   MPI_DOUBLE_INT, mpi_env.comm());
	} else {
	    MPI_Allgatherv(0, 0, MPI_DOUBLE_INT,
			   rbuf.begin(), dbuf_sz.begin(), dbuf_disp.begin(),
			   MPI_DOUBLE_INT, mpi_env.comm());
	}

	store_buffer(data, rbuf, adj_matrix);
    }

    std::cout << master << "reading done: " << timer(TIMER_GET) << std::endl;

    return true;
} // rw_read_adj
