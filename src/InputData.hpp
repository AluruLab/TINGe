/***
 *  $Id: InputData.hpp 597 2010-06-23 01:06:56Z zola $
 **
 *  File: InputData.hpp
 *  Created: Feb 04, 2008
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

#ifndef INPUT_DATA_HPP
#define INPUT_DATA_HPP

#include <jaz/numeric_table.hpp>
#include <set>
#include <string>
#include <vector>


struct InputData {
    InputData() : boot(0), etime(0), g_nrow(0), g_ncol(0), first_row(0),
		  last_row(0), has_TF(false), has_null_dist(false), mi_mean(0.0),
		  mi_sd(0.0), nd_sum(0.0), nd_sum_sqr(0.0), nd_mean(0.0), nd_sd(0.0),
		  nd_S(0), r_lim(0), mi_thr_used(0.0), fdr(-1.0) { }


    typedef double value_type;


    unsigned int boot;
    unsigned int etime;

    unsigned int g_nrow;
    unsigned int g_ncol;

    unsigned int num_rows;
    unsigned int first_row;
    unsigned int last_row;

    // local expression data
    jaz::numeric_table<value_type> exp_table;

    // probe names
    std::vector<std::string> probes;

    // TFs names
    std::set<std::string> tfs;

    // local TFs
    bool has_TF;
    std::vector<bool> isTF;

    // rank transformed local data (row-wise)
    jaz::plain_array<unsigned int> rank_array;

    // MI stats
    double mi_mean;
    double mi_sd;

    // null distribution
    bool has_null_dist;
    jaz::plain_array<value_type> null_dist;

    // sum, sum of squares, mean and sd of null distribution
    double nd_sum;
    double nd_sum_sqr;
    double nd_mean;
    double nd_sd;

    // null distribution size
    unsigned int nd_S;
    unsigned int r_lim;

    // actual MI threshold used
    double mi_thr_used;

    // FDR, if less than 0.0 it has not been set
    double fdr;

}; // struct InputData

#endif // INPUT_DATA_HPP
