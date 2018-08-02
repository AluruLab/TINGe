/***
 *  $Id: read_tfs.cpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: read_tfs.cpp
 *  Created: Apr 23, 2008
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

#ifndef READ_TFS_CPP
#define READ_TFS_CPP

#include "read_tfs.hpp"


bool read_tfs(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
	      InputData& data) {
    TIMER_START;
    std::cout << master << "* reading TFs list..." << std::endl;

    MPI_File fh;
    MPI_Status mpi_stat;

    int res = MPI_File_open(mpi_env.comm(),
			    const_cast<char*>(app_conf.input_tf.c_str()),
			    MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    if (res != MPI_SUCCESS) {
	std::cout << master << error << "can't open file " << app_conf.input_tf
		  << std::endl;
	return false;
    }

    unsigned int buf_sz = jaz::file_size(app_conf.input_tf.c_str());

    if (buf_sz == 0) {
	MPI_File_close(&fh);
	return true;
    }

    char* buf = new char[buf_sz];

    MPI_File_seek(fh, 0, MPI_SEEK_SET);
    MPI_File_read(fh, buf, buf_sz, MPI_CHAR, &mpi_stat);

    MPI_File_close(&fh);

    std::string s;
    std::istringstream ss;

    if (app_conf.has_pubsetbuf == true) ss.rdbuf()->pubsetbuf(buf, buf_sz);
    else ss.str(std::string(buf, buf_sz));

    std::copy(jaz::getline_iterator<>(ss), jaz::getline_iterator<>(),
	      std::inserter(data.tfs, data.tfs.end()));

    delete[] buf;

    std::cout << master << "number of TFs from file: " << data.tfs.size() << '\n';
    std::cout << master << "reading done: " << timer(TIMER_GET) << std::endl;

    return true;
} // read_tfs

#endif // READ_TFS_CPP
