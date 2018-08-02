/***
 *  $Id: read_tfs.hpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: read_tfs.hpp
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

#ifndef READ_TFS_HPP
#define READ_TFS_HPP

#include <mpix/MPI_env.hpp>
#include <jaz/iterator_add.hpp>
#include <jaz/sys_tools.hpp>
#include "AppConfig.hpp"
#include "InputData.hpp"
#include "utility.hpp"
#include <sstream>


bool read_tfs(const mpix::MPI_env& mpi_env, const AppConfig& app_conf,
	      InputData& data);

#endif // READ_TFS_HPP
